////////////////////////////////////
#ifndef KINEMATIC_H
#define KINEMATIC_H

////////////////////////////////////
#include "star.h"
#include "centroid.h"

#include <iostream>
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

        double minX = -8.0;             // kpc
        double maxX =  8.0;             // kpc
        double minY = -8.0;             // kpc
        double maxY =  8.0;             // kpc
        double minZ =  0.0;             // kpc
        double maxZ =  0.0;             // kpc
        double step = 0.25;             // kpc

        double starsRegionRadius = 1.0; // kpc
    };

    ////////////////////////////////

    std::vector<Centroid> calcCentroids(
            const std::list<Star>& allStars,
            const ConfigProcessing& config);
}

////////////////////////////////////
#endif // KINEMATIC_H
