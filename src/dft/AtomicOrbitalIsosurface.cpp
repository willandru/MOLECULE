#include "AtomicOrbitalIsosurface.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace
{

constexpr double EPSILON = 1.0e-12;

struct Tetrahedron
{
    int a;
    int b;
    int c;
    int d;
};

constexpr Tetrahedron TETRAHEDRA[6] =
{
    {0, 5, 1, 6},
    {0, 1, 2, 6},
    {0, 2, 3, 6},
    {0, 3, 7, 6},
    {0, 7, 4, 6},
    {0, 4, 5, 6}
};

}

AtomicOrbitalIsosurface::Result
AtomicOrbitalIsosurface::generate(
    const AtomicOrbital3D::Grid& grid,
    double isovalue
) const
{
    if (grid.nx < 2 ||
        grid.ny < 2 ||
        grid.nz < 2) {

        throw std::invalid_argument(
            "La malla 3D debe tener al menos 2 puntos por eje."
        );
    }

    const std::size_t expectedSize =
        grid.nx *
        grid.ny *
        grid.nz;

    if (grid.values.size() != expectedSize) {
        throw std::invalid_argument(
            "El tamano de los valores no coincide con la malla 3D."
        );
    }

    if (isovalue <= 0.0) {
        throw std::invalid_argument(
            "El isovalor debe ser mayor que cero."
        );
    }

    Result result;

    result.positive =
        generateSurface(
            grid,
            isovalue
        );

    result.negative =
        generateSurface(
            grid,
            -isovalue
        );

    return result;
}

AtomicOrbitalIsosurface::Surface
AtomicOrbitalIsosurface::generateSurface(
    const AtomicOrbital3D::Grid& grid,
    double isovalue
) const
{
    Surface surface;

    for (std::size_t z = 0;
         z < grid.nz - 1;
         ++z) {

        for (std::size_t y = 0;
             y < grid.ny - 1;
             ++y) {

            for (std::size_t x = 0;
                 x < grid.nx - 1;
                 ++x) {

                Point cube[8];

                const std::size_t x0 = x;
                const std::size_t x1 = x + 1;

                const std::size_t y0 = y;
                const std::size_t y1 = y + 1;

                const std::size_t z0 = z;
                const std::size_t z1 = z + 1;

                const double dx =
                    (grid.xmax - grid.xmin) /
                    static_cast<double>(grid.nx - 1);

                const double dy =
                    (grid.ymax - grid.ymin) /
                    static_cast<double>(grid.ny - 1);

                const double dz =
                    (grid.zmax - grid.zmin) /
                    static_cast<double>(grid.nz - 1);

                cube[0] = {
                    grid.xmin + x0 * dx,
                    grid.ymin + y0 * dy,
                    grid.zmin + z0 * dz,
                    grid.values[
                        gridIndex(
                            x0,
                            y0,
                            z0,
                            grid.nx,
                            grid.ny
                        )
                    ]
                };

                cube[1] = {
                    grid.xmin + x1 * dx,
                    grid.ymin + y0 * dy,
                    grid.zmin + z0 * dz,
                    grid.values[
                        gridIndex(
                            x1,
                            y0,
                            z0,
                            grid.nx,
                            grid.ny
                        )
                    ]
                };

                cube[2] = {
                    grid.xmin + x1 * dx,
                    grid.ymin + y1 * dy,
                    grid.zmin + z0 * dz,
                    grid.values[
                        gridIndex(
                            x1,
                            y1,
                            z0,
                            grid.nx,
                            grid.ny
                        )
                    ]
                };

                cube[3] = {
                    grid.xmin + x0 * dx,
                    grid.ymin + y1 * dy,
                    grid.zmin + z0 * dz,
                    grid.values[
                        gridIndex(
                            x0,
                            y1,
                            z0,
                            grid.nx,
                            grid.ny
                        )
                    ]
                };

                cube[4] = {
                    grid.xmin + x0 * dx,
                    grid.ymin + y0 * dy,
                    grid.zmin + z1 * dz,
                    grid.values[
                        gridIndex(
                            x0,
                            y0,
                            z1,
                            grid.nx,
                            grid.ny
                        )
                    ]
                };

                cube[5] = {
                    grid.xmin + x1 * dx,
                    grid.ymin + y0 * dy,
                    grid.zmin + z1 * dz,
                    grid.values[
                        gridIndex(
                            x1,
                            y0,
                            z1,
                            grid.nx,
                            grid.ny
                        )
                    ]
                };

                cube[6] = {
                    grid.xmin + x1 * dx,
                    grid.ymin + y1 * dy,
                    grid.zmin + z1 * dz,
                    grid.values[
                        gridIndex(
                            x1,
                            y1,
                            z1,
                            grid.nx,
                            grid.ny
                        )
                    ]
                };

                cube[7] = {
                    grid.xmin + x0 * dx,
                    grid.ymin + y1 * dy,
                    grid.zmin + z1 * dz,
                    grid.values[
                        gridIndex(
                            x0,
                            y1,
                            z1,
                            grid.nx,
                            grid.ny
                        )
                    ]
                };

                polygonizeCell(
                    cube,
                    isovalue,
                    surface
                );
            }
        }
    }

    calculateNormals(surface);

    return surface;
}

AtomicOrbitalIsosurface::EdgeVertex
AtomicOrbitalIsosurface::interpolate(
    const Point& a,
    const Point& b,
    double isovalue
) const
{
    if (std::abs(isovalue - a.value) < EPSILON) {
        return {
            a.x,
            a.y,
            a.z
        };
    }

    if (std::abs(isovalue - b.value) < EPSILON) {
        return {
            b.x,
            b.y,
            b.z
        };
    }

    const double denominator =
        b.value - a.value;

    if (std::abs(denominator) < EPSILON) {
        return {
            0.5 * (a.x + b.x),
            0.5 * (a.y + b.y),
            0.5 * (a.z + b.z)
        };
    }

    double t =
        (isovalue - a.value) /
        denominator;

    t = std::clamp(
        t,
        0.0,
        1.0
    );

    return {
        a.x + t * (b.x - a.x),
        a.y + t * (b.y - a.y),
        a.z + t * (b.z - a.z)
    };
}

void AtomicOrbitalIsosurface::polygonizeCell(
    const Point cube[8],
    double isovalue,
    Surface& surface
) const
{
    for (const Tetrahedron& tetra : TETRAHEDRA) {

        const Point points[4] =
        {
            cube[tetra.a],
            cube[tetra.b],
            cube[tetra.c],
            cube[tetra.d]
        };

        bool inside[4];

        int insideCount = 0;

        for (int i = 0; i < 4; ++i) {

            inside[i] =
                points[i].value >= isovalue;

            if (inside[i]) {
                ++insideCount;
            }
        }

        if (insideCount == 0 ||
            insideCount == 4) {

            continue;
        }

        std::vector<EdgeVertex> intersections;

        const int edges[6][2] =
        {
            {0, 1},
            {1, 2},
            {2, 3},
            {3, 0},
            {0, 2},
            {1, 3}
        };

        for (const auto& edge : edges) {

            const int a = edge[0];
            const int b = edge[1];

            if (inside[a] == inside[b]) {
                continue;
            }

            intersections.push_back(
                interpolate(
                    points[a],
                    points[b],
                    isovalue
                )
            );
        }

        if (intersections.size() == 3) {

            const unsigned int base =
                static_cast<unsigned int>(
                    surface.vertices.size()
                );

            for (const EdgeVertex& vertex :
                 intersections) {

                surface.vertices.push_back({
                    vertex.x,
                    vertex.y,
                    vertex.z,
                    0.0,
                    0.0,
                    0.0
                });
            }

            surface.indices.push_back(base);
            surface.indices.push_back(base + 1);
            surface.indices.push_back(base + 2);
        }
        else if (intersections.size() == 4) {

            const unsigned int base =
                static_cast<unsigned int>(
                    surface.vertices.size()
                );

            for (const EdgeVertex& vertex :
                 intersections) {

                surface.vertices.push_back({
                    vertex.x,
                    vertex.y,
                    vertex.z,
                    0.0,
                    0.0,
                    0.0
                });
            }

            surface.indices.push_back(base);
            surface.indices.push_back(base + 1);
            surface.indices.push_back(base + 2);

            surface.indices.push_back(base);
            surface.indices.push_back(base + 2);
            surface.indices.push_back(base + 3);
        }
    }
}

std::size_t AtomicOrbitalIsosurface::gridIndex(
    std::size_t x,
    std::size_t y,
    std::size_t z,
    std::size_t nx,
    std::size_t ny
) const
{
    return
        x +
        nx * (
            y +
            ny * z
        );
}

void AtomicOrbitalIsosurface::calculateNormals(
    Surface& surface
) const
{
    for (Vertex& vertex : surface.vertices) {

        vertex.nx = 0.0;
        vertex.ny = 0.0;
        vertex.nz = 0.0;
    }

    for (std::size_t i = 0;
         i + 2 < surface.indices.size();
         i += 3) {

        const unsigned int ia =
            surface.indices[i];

        const unsigned int ib =
            surface.indices[i + 1];

        const unsigned int ic =
            surface.indices[i + 2];

        Vertex& a =
            surface.vertices[ia];

        Vertex& b =
            surface.vertices[ib];

        Vertex& c =
            surface.vertices[ic];

        const double ux =
            b.x - a.x;

        const double uy =
            b.y - a.y;

        const double uz =
            b.z - a.z;

        const double vx =
            c.x - a.x;

        const double vy =
            c.y - a.y;

        const double vz =
            c.z - a.z;

        const double nx =
            uy * vz -
            uz * vy;

        const double ny =
            uz * vx -
            ux * vz;

        const double nz =
            ux * vy -
            uy * vx;

        a.nx += nx;
        a.ny += ny;
        a.nz += nz;

        b.nx += nx;
        b.ny += ny;
        b.nz += nz;

        c.nx += nx;
        c.ny += ny;
        c.nz += nz;
    }

    for (Vertex& vertex :
         surface.vertices) {

        const double length =
            std::sqrt(
                vertex.nx * vertex.nx +
                vertex.ny * vertex.ny +
                vertex.nz * vertex.nz
            );

        if (length > EPSILON) {

            vertex.nx /= length;
            vertex.ny /= length;
            vertex.nz /= length;
        }
    }
}