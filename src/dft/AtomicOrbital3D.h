#pragma once

#include "AtomicOrbitalAngular.h"
#include "DFTData.h"
#include "RadialGrid.h"

#include <cstddef>
#include <vector>

class AtomicOrbital3D
{
public:

    struct Grid
    {
        std::size_t nx = 0;
        std::size_t ny = 0;
        std::size_t nz = 0;

        double xmin = 0.0;
        double xmax = 0.0;

        double ymin = 0.0;
        double ymax = 0.0;

        double zmin = 0.0;
        double zmax = 0.0;

        std::vector<double> values;
    };

public:

    AtomicOrbital3D() = default;

    double evaluate(
        const RadialGrid& radialGrid,
        const AtomicOrbital& orbital,
        AtomicOrbitalAngularType angularType,
        double x,
        double y,
        double z
    ) const;

    Grid sample(
        const RadialGrid& radialGrid,
        const AtomicOrbital& orbital,
        AtomicOrbitalAngularType angularType,
        std::size_t pointsPerAxis,
        double extent
    ) const;

    Grid sample(
        const RadialGrid& radialGrid,
        const AtomicOrbital& orbital,
        int orbitalIndex,
        std::size_t pointsPerAxis,
        double extent
    ) const;

private:

    double radialValue(
        const RadialGrid& radialGrid,
        const AtomicOrbital& orbital,
        double radius
    ) const;

    std::size_t index(
        std::size_t ix,
        std::size_t iy,
        std::size_t iz,
        std::size_t nx,
        std::size_t ny
    ) const;
};