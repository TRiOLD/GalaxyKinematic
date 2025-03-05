////////////////////////////////////
#include "calccatalog.h"

#include <cmath>
#include <vector>
//#include <thread>

#include <TException.h>
#include <TFiletable.h>
#include <TLogger.h>
using namespace TRiOLD;

////////////////////////////////////
double _d2r(double angle)
{
    return angle / 180.0 * M_PI;
}

std::list<Star> CalcCatalog::readCatalog(const std::string& filepath,
                                         const ConfigTable& config)
{
    Filetable file;
    switch (config.fileType)
    {
        case 1: file.open(filepath, Filetable::IN, Filetable::TSV); break;
        case 2: file.open(filepath, Filetable::IN, Filetable::CSV); break;
        default: file.open(filepath);
    }

    if(!file.isOpen())
        Exception("Catalog::readCatalog() !file.isOpen()");

    Filetable::Row row;
    if(config.isWithHeader)
        file >> row;

    std::list<Star> res;
    while(file >> row)
    {
        Star::Equatorial EqC(
            (config.coordInRad) ? row.at(config.cn_ra).toDouble() :
                                  _d2r(row.at(config.cn_ra).toDouble()),
            (config.coordInRad) ? row.at(config.cn_de).toDouble() :
                                  _d2r(row.at(config.cn_de).toDouble()),
            1.0 / row.at(config.cn_plx).toDouble());
        Star::Equatorial EqV(
            row.at(config.cn_muRa).toDouble(),
            row.at(config.cn_muDe).toDouble(),
            row.at(config.cn_Vr).toDouble());
        res.push_back({EqC, EqV});
    }
    file.close();
    return res;
}

std::list<Star> CalcCatalog::readCatalog_agreedStruct(const std::string& filepath)
{
    Filetable file(filepath);
    if(!file.isOpen())
        Exception("Catalog::readCatalog() !file.isOpen()");

    Filetable::Row row;
    file >> row;

    std::list<Star> res;
    while(file >> row)
    {
        Star::Cartesian GCC(
            row.at(0).toDouble(), row.at(1).toDouble(), row.at(2).toDouble());
        Star::Cartesian GCV(
            row.at(3).toDouble(), row.at(4).toDouble(), row.at(5).toDouble());
        res.push_back({GCC, GCV});
    }
    file.close();
    return res;
}

std::list<Star> CalcCatalog::createPixCatalog_agreedStruct(const std::list<Star> &stars,
                                                           const ConfigPixelization& config)
{
    unsigned int areaSizeX = (config.maxX - config.minX) / (config.pixelHalfWidth * 2.0) + 1;
    unsigned int areaSizeY = (config.maxY - config.minY) / (config.pixelHalfWidth * 2.0) + 1;
    unsigned int areaSizeZ = (config.maxZ - config.minZ) / (config.pixelHalfWidth * 2.0) + 1;

    typedef std::list<Star::Cartesian> PixelV;
    typedef std::vector<PixelV> PixelsVoX;
    typedef std::vector<std::vector<PixelV>> PixelsVXoY;
    typedef std::vector<std::vector<std::vector<PixelV>>> PixelsV3D;

    PixelsV3D pixelsV(areaSizeZ, PixelsVXoY(areaSizeY, PixelsVoX(areaSizeX)));
    for(const Star &s : stars) {
        Star::Cartesian GCC = s.getGCC();
        if(GCC.x > config.maxX + config.pixelHalfWidth ||
           GCC.x < config.minX - config.pixelHalfWidth ||
           GCC.y > config.maxY + config.pixelHalfWidth ||
           GCC.y < config.minY - config.pixelHalfWidth ||
           GCC.z > config.maxZ + config.pixelHalfWidth ||
           GCC.z < config.minZ - config.pixelHalfWidth) {
            continue;
        }
        int i = (GCC.x - config.minX - config.pixelHalfWidth) / (config.pixelHalfWidth * 2.0);
        int j = (GCC.y - config.minY - config.pixelHalfWidth) / (config.pixelHalfWidth * 2.0);
        int k = (GCC.z - config.minZ - config.pixelHalfWidth) / (config.pixelHalfWidth * 2.0);
        pixelsV.at(k).at(j).at(i).push_back(s.getGCV());
    }

    std::list<Star> res;
    for(unsigned int k = 0; k < areaSizeZ; ++k) {
        for(unsigned int j = 0; j < areaSizeY; ++j) {
            for(unsigned int i = 0; i < areaSizeX; ++i) {
                const PixelV &pixelV = pixelsV.at(k).at(j).at(i);
                if(pixelV.size() >= config.minStarCount) {
                    Star::Cartesian averV;
                    for(const Star::Cartesian &subv : pixelV) {
                        averV.x += subv.x / (double)pixelV.size();
                        averV.y += subv.y / (double)pixelV.size();
                        averV.z += subv.z / (double)pixelV.size();
                    }
                    double X = (double)i * config.pixelHalfWidth * 2.0 + config.minX;
                    double Y = (double)j * config.pixelHalfWidth * 2.0 + config.minY;
                    double Z = (double)k * config.pixelHalfWidth * 2.0 + config.minZ;
                    res.push_back(Star(Star::Cartesian(X, Y, Z), averV));
                }
            }
        }
    }
    return res;
}

void CalcCatalog::writeCatalog_agreedStruct(const std::string& filepath,
                                            const std::list<Star>& stars)
{
    Filetable file(filepath, Filetable::OUT, Filetable::CSV);
    if(!file.isOpen())
        Exception("Catalog::readCatalog() !file.isOpen()");

    file << Filetable::Row({
        {"X[kpc]"}, {"Y[kpc]"}, {"Z[kpc]"},
        {"VX[km/s]"}, {"VY[km/s]"}, {"VZ[km/s]"}});

    for(const auto& it : stars)
    {
        Star::Cartesian GCC = it.getGCC();
        Star::Cartesian GCV = it.getGCV();
        file << Filetable::Row(
            {{GCC.x, 12}, {GCC.y, 12}, {GCC.z, 12},
             {GCV.x, 12}, {GCV.y, 12}, {GCV.z, 12}});
    }
    file.close();
}

////////////////////////////////////

void CalcCatalog::writeCentroids(const std::string& filepath,
                                 const std::vector<Centroid>& centroids,
                                 const ConfigConstants& config)
{
    Filetable file(filepath, Filetable::OUT, Filetable::CSV);
    if(!file.isOpen())
        Exception("Catalog::readCatalog() !file.isOpen()");

    file << Filetable::Row({
        {"starsAmount"},
        {"X[kpc]"}, {"Y[kpc]"}, {"Z[kpc]"},
        {"VX[km/s]"}, {"VY[km/s]"}, {"VZ[km/s]"},
        {"R[kpc]"}, {"theta[rad]"}, {"Z[kpc]"},
        {"VR[km/s]"}, {"Vtheta[km/s]"}, {"VZ[km/s]"},

        {"u[km/s]"}, {"v[km/s]"}, {"w[km/s]"},
        {"w1[km/s/kpc]"}, {"w2[km/s/kpc]"}, {"w3[km/s/kpc]"},
        {"mp12[km/s/kpc]"}, {"mp23[km/s/kpc]"}, {"mp13[km/s/kpc]"},
        {"mp11[km/s/kpc]"}, {"mp22[km/s/kpc]"}, {"mp33[km/s/kpc]"},

        {"u_err[km/s]"}, {"v_err[km/s]"}, {"w_err[km/s]"},
        {"w1_err[km/s/kpc]"}, {"w2_err[km/s/kpc]"}, {"w3_err[km/s/kpc]"},
        {"mp12_err[km/s/kpc]"}, {"mp23_err[km/s/kpc]"}, {"mp13_err[km/s/kpc]"},
        {"mp11_err[km/s/kpc]"}, {"mp22_err[km/s/kpc]"}, {"mp33_err[km/s/kpc]"},

        {"u_rot[km/s]"}, {"v_rot[km/s]"}, {"w_rot[km/s]"},
        {"w1_rot[km/s/kpc]"}, {"w2_rot[km/s/kpc]"}, {"w3_rot[km/s/kpc]"},
        {"mp12_rot[km/s/kpc]"}, {"mp23_rot[km/s/kpc]"}, {"mp13_rot[km/s/kpc]"},
        {"mp11_rot[km/s/kpc]"}, {"mp22_rot[km/s/kpc]"}, {"mp33_rot[km/s/kpc]"}});

    for(const auto& it : centroids)
    {
        Centroid::Cartesian GCC = it.getGCC();
        Centroid::Cartesian GCV = it.getGCV();
        Centroid::Cylindrical GCCC, GCCV;
        it.getcalcGCCCandV(GCCC, GCCV,
            config.R_Sun, {config.VX_Sun, config.VY_Sun, config.VZ_Sun});
        Centroid::KinematicParameters KPs = it.getKPs();
        Centroid::KinematicParameters KPsErr = it.getKPsErr();
        Centroid::KinematicParameters KPsRot = it.getcalcKPs_localRot(config.R_Sun);

        const int p = 12; // precision
        file << Filetable::Row({
            {(int)it.getStarsAmount()},
            {GCC.x, p}, {GCC.y, p}, {GCC.z, p},
            {GCV.x, p}, {GCV.y, p}, {GCV.z, p},
            {GCCC.R, p}, {GCCC.theta, p}, {GCCC.Z, p},
            {GCCV.R, p}, {GCCV.theta, p}, {GCCV.Z, p},

            {KPs.u, p}, {KPs.v, p}, {KPs.w, p},
            {KPs.w1, p}, {KPs.w2, p}, {KPs.w3, p},
            {KPs.mp12, p}, {KPs.mp23, p}, {KPs.mp13, p},
            {KPs.mp11, p}, {KPs.mp22, p}, {KPs.mp33, p},

            {KPsErr.u, p}, {KPsErr.v, p}, {KPsErr.w, p},
            {KPsErr.w1, p}, {KPsErr.w2, p}, {KPsErr.w3, p},
            {KPsErr.mp12, p}, {KPsErr.mp23, p}, {KPsErr.mp13, p},
            {KPsErr.mp11, p}, {KPsErr.mp22, p}, {KPsErr.mp33, p},

            {KPsRot.u, p}, {KPsRot.v, p}, {KPsRot.w, p},
            {KPsRot.w1, p}, {KPsRot.w2, p}, {KPsRot.w3, p},
            {KPsRot.mp12, p}, {KPsRot.mp23, p}, {KPsRot.mp13, p},
            {KPsRot.mp11, p}, {KPsRot.mp22, p}, {KPsRot.mp33, p}});
    }
    file.close();
}
