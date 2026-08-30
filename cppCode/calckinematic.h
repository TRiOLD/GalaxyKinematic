////////////////////////////////////
#ifndef KINEMATIC_H
#define KINEMATIC_H

////////////////////////////////////
#include "star.h"
#include "centroid.h"

#include <vector>
#include <list>

////////////////////////////////////
namespace CalcKinematic
{
    ////////////////////////////////
    struct ConfigProcessing
    {
        unsigned int threadsAmount = 16;
        unsigned int RAMlimit = 16000;  // MB

        int gridType = 0; // 0 = cartesian, 1 = cylindrical

        // cartesian grid params
        double minX = -8.0;             // kpc
        double maxX =  8.0;             // kpc
        double stepX = 0.1;             // kpc
        double minY = -8.0;             // kpc
        double maxY =  8.0;             // kpc
        double stepY = 0.1;             // kpc

        // cylindrical grid params
        double minR =  0.0;             // kpc
        double maxR = 25.0;             // kpc
        double stepR = 0.1;             // kpc
        double shiftR = 8.28;           // kpc
        double minTheta =  90.0;        // deg
        double maxTheta = 270.0;        // deg
        double stepTheta = 1.0;         // deg

        // common grid options
        double minZ =  0.0;             // kpc
        double maxZ =  0.0;             // kpc
        double stepZ = 0.5;             // kpc

        double starsRegionRadius = 1.0; // kpc
    };

    ////////////////////////////////

    std::vector<Centroid> calcCentroids(
            const std::list<Star> &allStars,
            const ConfigProcessing &config);
}

////////////////////////////////////
#endif // KINEMATIC_H
