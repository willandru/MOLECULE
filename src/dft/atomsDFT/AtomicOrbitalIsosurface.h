#pragma once

#include "AtomicOrbital3D.h"

#include <cstddef>
#include <vector>

class AtomicOrbitalIsosurface
{
public:

    struct Vertex
    {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;

        double nx = 0.0;
        double ny = 0.0;
        double nz = 0.0;
    };

    struct Surface
    {
        std::vector<Vertex> vertices;
        std::vector<unsigned int> indices;
    };

    struct Result
    {
        Surface positive;
        Surface negative;
    };

public:

    AtomicOrbitalIsosurface() = default;

    Result generate(
        const AtomicOrbital3D::Grid& grid,
        double isovalue
    ) const;

private:

    struct Point
    {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;

        double value = 0.0;
    };

    struct EdgeVertex
    {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
    };

private:

    Surface generateSurface(
        const AtomicOrbital3D::Grid& grid,
        double isovalue
    ) const;

    EdgeVertex interpolate(
        const Point& a,
        const Point& b,
        double isovalue
    ) const;

    void polygonizeCell(
        const Point cube[8],
        double isovalue,
        Surface& surface
    ) const;

    std::size_t gridIndex(
        std::size_t x,
        std::size_t y,
        std::size_t z,
        std::size_t nx,
        std::size_t ny
    ) const;

    void calculateNormals(
        Surface& surface
    ) const;
};