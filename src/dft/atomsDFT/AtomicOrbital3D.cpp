#include "AtomicOrbital3D.h"

#include <cmath>
#include <stdexcept>
#include <vector>

namespace
{

double interpolateRadialFunction(
    const RadialGrid& radialGrid,
    const AtomicOrbital& orbital,
    double radius
)
{
    const std::vector<double>& coordinates =
        radialGrid.coordinates();

    const std::vector<double>& radial =
        orbital.u;

    if (coordinates.empty() || radial.empty()) {
        return 0.0;
    }

    if (coordinates.size() != radial.size()) {
        throw std::invalid_argument(
            "La malla radial y el orbital tienen tamanos incompatibles."
        );
    }

    if (radius <= 0.0) {
        return 0.0;
    }

    if (radius <= coordinates.front()) {

        const double r =
            coordinates.front();

        if (r <= 0.0) {
            return 0.0;
        }

        return radial.front() / r;
    }

    if (radius >= coordinates.back()) {
        return 0.0;
    }

    const double spacing =
        radialGrid.spacing();

    if (spacing <= 0.0) {
        return 0.0;
    }

    const double position =
        radius / spacing - 1.0;

    const std::size_t lowerIndex =
        static_cast<std::size_t>(
            std::floor(position)
        );

    const std::size_t upperIndex =
        lowerIndex + 1;

    if (upperIndex >= radial.size()) {
        return 0.0;
    }

    const double r0 =
        coordinates[lowerIndex];

    const double r1 =
        coordinates[upperIndex];

    const double u0 =
        radial[lowerIndex];

    const double u1 =
        radial[upperIndex];

    const double interval =
        r1 - r0;

    if (interval <= 0.0) {
        return u0 / r0;
    }

    const double fraction =
        (radius - r0) / interval;

    const double u =
        u0 +
        fraction * (u1 - u0);

    return u / radius;
}

}

double AtomicOrbital3D::evaluate(
    const RadialGrid& radialGrid,
    const AtomicOrbital& orbital,
    AtomicOrbitalAngularType angularType,
    double x,
    double y,
    double z
) const
{
    const double radius =
        std::sqrt(
            x * x +
            y * y +
            z * z
        );

    if (radius <= 0.0) {

        if (angularType ==
            AtomicOrbitalAngularType::S) {

            const double firstRadius =
                radialGrid.coordinates().front();

            return
                radialValue(
                    radialGrid,
                    orbital,
                    firstRadius
                ) *
                AtomicOrbitalAngular::evaluateNormalized(
                    angularType,
                    0.0,
                    0.0,
                    0.0
                );
        }

        return 0.0;
    }

    const double radial =
        radialValue(
            radialGrid,
            orbital,
            radius
        );

    if (radial == 0.0) {
        return 0.0;
    }

    const double angular =
        AtomicOrbitalAngular::evaluateNormalized(
            angularType,
            x,
            y,
            z
        );

    return radial * angular;
}

AtomicOrbital3D::Grid AtomicOrbital3D::sample(
    const RadialGrid& radialGrid,
    const AtomicOrbital& orbital,
    AtomicOrbitalAngularType angularType,
    std::size_t pointsPerAxis,
    double extent
) const
{
    if (pointsPerAxis < 2) {
        throw std::invalid_argument(
            "pointsPerAxis debe ser al menos 2."
        );
    }

    if (extent <= 0.0) {
        throw std::invalid_argument(
            "extent debe ser mayor que cero."
        );
    }

    Grid grid;

    grid.nx = pointsPerAxis;
    grid.ny = pointsPerAxis;
    grid.nz = pointsPerAxis;

    grid.xmin = -extent;
    grid.xmax = extent;

    grid.ymin = -extent;
    grid.ymax = extent;

    grid.zmin = -extent;
    grid.zmax = extent;

    grid.values.resize(
        pointsPerAxis *
        pointsPerAxis *
        pointsPerAxis,
        0.0
    );

    const double denominator =
        static_cast<double>(
            pointsPerAxis - 1
        );

    const double dx =
        (grid.xmax - grid.xmin) /
        denominator;

    const double dy =
        (grid.ymax - grid.ymin) /
        denominator;

    const double dz =
        (grid.zmax - grid.zmin) /
        denominator;

    for (std::size_t iz = 0;
         iz < pointsPerAxis;
         ++iz) {

        const double z =
            grid.zmin +
            static_cast<double>(iz) * dz;

        for (std::size_t iy = 0;
             iy < pointsPerAxis;
             ++iy) {

            const double y =
                grid.ymin +
                static_cast<double>(iy) * dy;

            for (std::size_t ix = 0;
                 ix < pointsPerAxis;
                 ++ix) {

                const double x =
                    grid.xmin +
                    static_cast<double>(ix) * dx;

                grid.values[
                    index(
                        ix,
                        iy,
                        iz,
                        grid.nx,
                        grid.ny
                    )
                ] =
                    evaluate(
                        radialGrid,
                        orbital,
                        angularType,
                        x,
                        y,
                        z
                    );
            }
        }
    }

    return grid;
}

AtomicOrbital3D::Grid AtomicOrbital3D::sample(
    const RadialGrid& radialGrid,
    const AtomicOrbital& orbital,
    int orbitalIndex,
    std::size_t pointsPerAxis,
    double extent
) const
{
    const AtomicOrbitalAngularType angularType =
        AtomicOrbitalAngular::typeFromQuantumNumbers(
            orbital.l,
            orbitalIndex
        );

    return sample(
        radialGrid,
        orbital,
        angularType,
        pointsPerAxis,
        extent
    );
}

double AtomicOrbital3D::radialValue(
    const RadialGrid& radialGrid,
    const AtomicOrbital& orbital,
    double radius
) const
{
    return interpolateRadialFunction(
        radialGrid,
        orbital,
        radius
    );
}

std::size_t AtomicOrbital3D::index(
    std::size_t ix,
    std::size_t iy,
    std::size_t iz,
    std::size_t nx,
    std::size_t ny
) const
{
    return
        ix +
        nx * (
            iy +
            ny * iz
        );
}