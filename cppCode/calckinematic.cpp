////////////////////////////////////
#include "calckinematic.h"

#include <condition_variable>
#include <functional>
#include <atomic>
#include <mutex>
#include <queue>
#include <thread>

//ALGLIB Lib
#include <interpolation.h>
using namespace alglib;

#include <TException.h>
#include <TLogger.h>
using namespace TRiOLD;

////////////////////////////////////
namespace {
struct CentroidProcessingData {
    Centroid::Cartesian centroidGCC;
    std::size_t localStarsAmount;
    std::size_t estimatedRAM;
};

class ThreadPool
{
public:
    ThreadPool(size_t threads)
    {
        for (size_t i = 0; i < threads; ++i) {
            workers.emplace_back([this]() {
                while (true) {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lock(queue_mutex);
                        condition.wait(lock, [this]() {
                            return stop || !tasks.empty();
                        });
                        if (stop && tasks.empty()) {
                            return;
                        }
                        task = std::move(tasks.front());
                        tasks.pop();
                    }
                    task();
                    tasks_in_progress--;
                    condition.notify_all();
                }
            });
        }
    }

    void enqueue(std::function<void()> task)
    {
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            tasks.push(task);
            tasks_in_progress++;
        }
        condition.notify_one();
    }

    void wait()
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        condition.wait(lock, [this]() {
            return tasks.empty() && tasks_in_progress == 0;
        });
    }

    ~ThreadPool()
    {
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            stop = true;
        }
        condition.notify_all();

        for (std::thread &w : workers) {
            w.join();
        }
    }

private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;

    std::mutex queue_mutex;
    std::condition_variable condition;
    std::atomic<size_t> tasks_in_progress{0};
    bool stop = false;
};

////////////////////////////////////
class MemoryLimiter
{
public:
    explicit MemoryLimiter(size_t limitMB)
        : limit(limitMB), used(0)
    {}

    void acquire(size_t amount)
    {
        std::unique_lock lock(mtx);
        cv.wait(lock, [&]{
            return used + amount <= limit;
        });
        used += amount;
    }

    void release(size_t amount)
    {
        {
            std::lock_guard lock(mtx);
            used -= amount;
        }
        cv.notify_all();
    }

private:
    size_t limit;
    size_t used;

    std::mutex mtx;
    std::condition_variable cv;
};


////////////////////////////////////
bool _isInSphere(double radius, const Star::Cartesian &localGC)
{
    return std::sqrt(
        localGC.x*localGC.x +
        localGC.y*localGC.y +
        localGC.z*localGC.z) <= radius;
}

void _calcStarsAmount(
        std::size_t &starsAmount,
        const std::list<Star> &allStars,
        const Centroid::Cartesian &centroidGCC,
        double starsRegionRadius)
{
    starsAmount = 0;
    for (const Star &s : allStars) {
        if (_isInSphere(starsRegionRadius, s.getcalcGCC_local(centroidGCC))) {
            ++starsAmount;
        }
    }
}

void _selectStars(
        std::list<const Star *> &localStars_ptrs,
        const std::list<Star> &allStars,
        const Centroid::Cartesian &centroidGCC,
        double starsRegionRadius)
{
    localStars_ptrs.clear();
    for (const Star &s : allStars) {
        if (_isInSphere(starsRegionRadius, s.getcalcGCC_local(centroidGCC))) {
            localStars_ptrs.push_back(&s);
        }
    }
}

std::size_t _estimateThreadRAM(std::size_t localStarsSize)
{
    constexpr std::size_t bytesPerStar = 3072; // ~3 KB (with alglib matrix)
    constexpr std::size_t bytesPerMB = 1024 * 1024;
    return (localStarsSize * bytesPerStar) / bytesPerMB;
}

std::vector<CentroidProcessingData> _dataPreparation(
        const CalcKinematic::ConfigProcessing &config,
        const std::list<Star> &allStars)
{
    std::size_t cAmountX = (config.maxX - config.minX) / config.step + 1;
    std::size_t cAmountY = (config.maxY - config.minY) / config.step + 1;
    std::size_t cAmountZ = (config.maxZ - config.minZ) / config.step + 1;
    std::size_t centroidsAmount = cAmountX * cAmountY * cAmountZ;
    std::vector<CentroidProcessingData> res(centroidsAmount);

    std::size_t c = 0;
    std::size_t maxThreadsAmount = config.threadsAmount;
    if (maxThreadsAmount > centroidsAmount) {
        maxThreadsAmount = centroidsAmount;
    }
    ThreadPool pool(maxThreadsAmount);
    std::mutex maxEstimatedRAMMtx;
    size_t maxEstimatedRAM = 0;
    for (std::size_t i = 0; i < cAmountX; ++i)
        for (std::size_t j = 0; j < cAmountY; ++j)
            for (std::size_t k = 0; k < cAmountZ; ++k) {
                pool.enqueue([&, c, i, j, k]() {
                    res.at(c).centroidGCC = Centroid::Cartesian(
                        config.minX + config.step * i,
                        config.minY + config.step * j,
                        config.minZ + config.step * k);
                    _calcStarsAmount(res.at(c).localStarsAmount, allStars,
                                     res.at(c).centroidGCC, config.starsRegionRadius);
                    res.at(c).estimatedRAM = _estimateThreadRAM(res.at(c).localStarsAmount);
                    {
                        std::lock_guard lock(maxEstimatedRAMMtx);
                        if (maxEstimatedRAM < res.at(c).estimatedRAM) {
                            maxEstimatedRAM = res.at(c).estimatedRAM;
                        }
                    }
                });
                ++c;
            }
    pool.wait();
    if (maxEstimatedRAM > config.RAMlimit) {
        throw Exception("Needed minimum " + std::to_string(maxEstimatedRAM) + "MB RAM for processing");
    }
    LOG.writeInfo(std::to_string(centroidsAmount) + " centroids has been prepared.");
    return res;
}

void _createVariables( // by O-M kinematic model (GC)
        std::vector<double> &foos, std::vector<std::vector<double>> &polisVars,
        const Star::Cartesian &GCC, const Star::Cartesian &GCV)
{
    foos = {GCV.x, GCV.y, GCV.z};
    polisVars = {
        {-1.0, 0.0, 0.0, 0.0, GCC.z, -GCC.y, GCC.y, 0.0, GCC.z, GCC.x, 0.0, 0.0},
        {0.0, -1.0, 0.0, -GCC.z, 0.0, GCC.x, GCC.x, GCC.z, 0.0, 0.0, GCC.y, 0.0},
        {0.0, 0.0, -1.0, GCC.y, -GCC.x, 0.0, 0.0, GCC.y, GCC.x, 0.0, 0.0, GCC.z}
    };
}

Centroid _initNoCalcCentroid(const CentroidProcessingData& CPD)
{
    return Centroid(CPD.centroidGCC, {}, CPD.localStarsAmount, {}, Matrix<double>(3));
}

void _calcCentroid(
        Centroid &centroid,
        const CentroidProcessingData &CPD,
        double starsRegionRadius,
        const std::list<Star> &allStars)
{
    unsigned int starsAmount = CPD.localStarsAmount;
    if (starsAmount < 4) {
        centroid = _initNoCalcCentroid(CPD);
        return;
    }
    std::list<const Star *> localStars_ptrs;
    _selectStars(localStars_ptrs, allStars, CPD.centroidGCC, starsRegionRadius);
    real_1d_array y;
    real_2d_array fmatrix;
    y.setlength(3*starsAmount);
    fmatrix.setlength(3*starsAmount, 12);

    std::size_t e = 0;
    for (const Star *s : localStars_ptrs) {
        std::vector<double> foos;
        std::vector<std::vector<double>> polisVars;
        _createVariables(foos, polisVars,
            s->getcalcGCC_local(CPD.centroidGCC), s->getGCV());
        for (std::size_t r = 0; r < 3; ++r) {
            y[e+r] = foos[r];
            for (std::size_t c = 0; c < 12; ++c) {
                fmatrix[e+r][c] = polisVars[r][c];
            }
        }
        e += 3;
    }
    try {
        ae_int_t info;
        real_1d_array x;
        lsfitreport rep;
        lsfitlinear(y, fmatrix, info, x, rep);
        real_1d_array e = rep.errpar;

        centroid = Centroid(CPD.centroidGCC, {-x[0], -x[1], -x[2]}, starsAmount,
            {x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11]});
        centroid.setKPsErr(
            {e[0], e[1], e[2], e[3], e[4], e[5], e[6], e[7], e[8], e[9], e[10], e[11]});
    } catch (alglib::ap_error exc) {
        centroid = _initNoCalcCentroid(CPD);
        throw Exception(exc.msg, -30);
    }
}
}

std::vector<Centroid> CalcKinematic::calcCentroids(
        const std::list<Star> &allStars,
        const CalcKinematic::ConfigProcessing &config)
{
    std::vector<CentroidProcessingData> CPDs = _dataPreparation(config, allStars);
    std::size_t centroidsAmount = CPDs.size();
    std::vector<Centroid> res(centroidsAmount);

    std::size_t maxThreadsAmount = config.threadsAmount;
    if (maxThreadsAmount > centroidsAmount) {
        maxThreadsAmount = centroidsAmount;
    }
    ThreadPool pool(maxThreadsAmount);
    std::atomic<size_t> doneAnount{0};
    MemoryLimiter ramLimiter(config.RAMlimit);
    for (std::size_t c = 0; c < centroidsAmount; ++c) {
        pool.enqueue([&, c]() {
            ramLimiter.acquire(CPDs.at(c).estimatedRAM);
            try {
                _calcCentroid(res.at(c), CPDs.at(c), config.starsRegionRadius, allStars);
            } catch (...) {
                ramLimiter.release(CPDs.at(c).estimatedRAM);
                throw;
            }
            ramLimiter.release(CPDs.at(c).estimatedRAM);

            ++doneAnount;
            if (doneAnount % 100 == 0 || doneAnount == centroidsAmount) {
                LOG.writeInfo(std::to_string(doneAnount) + " / " + std::to_string(centroidsAmount) +
                    " centroids has been calculated.");
            }
        });
    }
    pool.wait();
    return res;
}

////////////////////////////////////
