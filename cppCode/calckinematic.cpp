////////////////////////////////////
#include "calckinematic.h"

#include <iostream>
#include <algorithm>
#include <functional>
#include <thread>

//ALGLIB Lib
#include <interpolation.h>
//#include <optimization.h>
//#include <specialfunctions.h>
//#include <statistics.h>
using namespace alglib;

#include <TException.h>
#include <TLogger.h>
using namespace TRiOLD;

////////////////////////////////////
struct CentroidProcessingData
{
    Centroid::Cartesian centroidGCC;
    std::size_t localStarsAmount;
};

////////////////////////////////////
bool _isInSphere(double radius, const Star::Cartesian& localGC)
{
    return std::sqrt(
        localGC.x*localGC.x +
        localGC.y*localGC.y +
        localGC.z*localGC.z) <= radius;
}

void _calcStarsAmount(
        std::size_t& starsAmount,
        const std::list<Star>& allStars,
        const Centroid::Cartesian &centroidGCC,
        double starsRegionRadius)
{
    starsAmount = 0;
    for(const auto& it : allStars)
        if(_isInSphere(starsRegionRadius, it.getcalcGCC_local(centroidGCC)))
            ++starsAmount;
}

void _selectStars(
        std::list<const Star*>& localStars_ptrs,
        const std::list<Star>& allStars,
        const Centroid::Cartesian &centroidGCC,
        double starsRegionRadius)
{
    localStars_ptrs.clear();
    for(const auto& it : allStars)
        if(_isInSphere(starsRegionRadius, it.getcalcGCC_local(centroidGCC)))
            localStars_ptrs.push_back(&it);
}

std::vector<CentroidProcessingData> _dataPreparation(
        const CalcKinematic::ConfigProcessing& config,
        const std::list<Star>& allStars)
{
    std::size_t cAmountX = (config.maxX - config.minX) / config.step + 1;
    std::size_t cAmountY = (config.maxY - config.minY) / config.step + 1;
    std::size_t cAmountZ = (config.maxZ - config.minZ) / config.step + 1;
    std::size_t cetroidsAmount = cAmountX * cAmountY * cAmountZ;
    std::vector<CentroidProcessingData> res(cetroidsAmount);

    std::size_t threadsAmount = config.threadsAmount;
    if(threadsAmount > cetroidsAmount) threadsAmount = cetroidsAmount;
    std::vector<std::thread> threads(threadsAmount);

    std::size_t c = 0, t = 0;
    for(std::size_t i = 0; i < cAmountX; ++i)
        for(std::size_t j = 0; j < cAmountY; ++j)
            for(std::size_t k = 0; k < cAmountZ; ++k)
            {
                res.at(c).centroidGCC = Centroid::Cartesian(
                    config.minX + config.step * i,
                    config.minY + config.step * j,
                    config.minZ + config.step * k);

                threads.at(t) = std::thread(
                    _calcStarsAmount, std::ref(res.at(c).localStarsAmount), std::cref(allStars),
                    std::cref(res.at(c).centroidGCC), config.starsRegionRadius);
                ++t;
                ++c;

                if(t == threadsAmount || c == cetroidsAmount)
                {
                    t = 0;
                    for(auto& th : threads)
                        if(th.joinable())
                            th.join();
                }
            }

    std::sort(res.begin(), res.end(),
        [](const CentroidProcessingData& CPD1, const CentroidProcessingData& CPD2)
        {return (CPD1.localStarsAmount > CPD2.localStarsAmount);});
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
        Centroid& centroid,
        const CentroidProcessingData& CPD,
        double starsRegionRadius,
        const std::list<Star>& allStars)
{
    unsigned int starsAmount = CPD.localStarsAmount;
    if(starsAmount < 4)
    {
        centroid = _initNoCalcCentroid(CPD);
        return;
    }

    std::list<const Star*> localStars_ptrs;
    _selectStars(localStars_ptrs, allStars,
                 CPD.centroidGCC, starsRegionRadius);

    real_1d_array y;
    real_2d_array fmatrix;
    y.setlength(3*starsAmount);
    fmatrix.setlength(3*starsAmount, 12);

    std::size_t e = 0;
    for(auto it : localStars_ptrs)
    {
        std::vector<double> foos;
        std::vector<std::vector<double>> polisVars;
        _createVariables(foos, polisVars,
            it->getcalcGCC_local(CPD.centroidGCC), it->getGCV());
        for(std::size_t r = 0; r < 3; ++r)
        {
            y[e+r] = foos[r];
            for(std::size_t c = 0; c < 12; ++c)
                fmatrix[e+r][c] = polisVars[r][c];
        }
        e += 3;
    }

    try
    {
        ae_int_t info;
        real_1d_array x;
        lsfitreport rep;
        lsfitlinear(y, fmatrix, info, x, rep);
        real_1d_array e = rep.errpar;

        centroid = Centroid(CPD.centroidGCC, {-x[0], -x[1], -x[2]}, starsAmount,
            {x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11]});
        centroid.setKPsErr(
            {e[0], e[1], e[2], e[3], e[4], e[5], e[6], e[7], e[8], e[9], e[10], e[11]});
    }
    catch (alglib::ap_error exc)
    {
        centroid = _initNoCalcCentroid(CPD);
        throw (Exception(exc.msg, -30));
    }
}

std::size_t _estimateStartRAM(std::size_t allStarsAmount)
{
    return std::size_t(allStarsAmount / 7500.0);
}

std::size_t _estimateThreadRAM(std::size_t localStarsSize)
{
    return std::size_t(localStarsSize / 250.0);
}

std::vector<Centroid> CalcKinematic::calcCentroids(
        const std::list<Star>& allStars,
        const CalcKinematic::ConfigProcessing& config)
{
    std::vector<CentroidProcessingData> CPDs = _dataPreparation(config, allStars);
    std::size_t cetroidsAmount = CPDs.size();
    std::vector<Centroid> res(cetroidsAmount);

    std::size_t maxThreadsAmount = config.threadsAmount;
    if(maxThreadsAmount > cetroidsAmount) maxThreadsAmount = cetroidsAmount;
    std::vector<std::thread> threads(maxThreadsAmount);

    std::size_t threadsRAM = 0;
    std::size_t startRAM = _estimateStartRAM(allStars.size());
    for(std::size_t c = 0, t = 0; c < cetroidsAmount;)
    {
        threadsRAM += _estimateThreadRAM(CPDs.at(c).localStarsAmount);
        bool isRAMoverflow = startRAM+threadsRAM >= config.RAMlimit;
        if(!isRAMoverflow && t < maxThreadsAmount)
        {
            threads.at(t) = std::thread(
                _calcCentroid, std::ref(res.at(c)), std::cref(CPDs.at(c)),
                config.starsRegionRadius, std::cref(allStars));
            ++t;
            ++c;
        }

        if(isRAMoverflow || t == 0)
            Exception("Kinematic::calcCentroids() Not enough RAM!");

        if(isRAMoverflow || t == maxThreadsAmount || c == cetroidsAmount)
        {
            t = 0;
            threadsRAM = 0;
            for(auto& th : threads)
                if(th.joinable())
                    th.join();
        }
    }
    return res;
}

////////////////////////////////////
