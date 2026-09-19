#include "MolecularOrbitalIsosurface.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace
{

constexpr double VALUE_EPSILON = 1.0e-14;

struct TetrahedronEdge
{
    int a;
    int b;
};

constexpr TetrahedronEdge TETRAHEDRON_EDGES[6] =
{
    { 0, 1 },
    { 1, 2 },
    { 2, 0 },
    { 0, 3 },
    { 1, 3 },
    { 2, 3 }
};

bool isInside(
    double value,
    double isovalue
)
{
    return value >= isovalue;
}

double interpolationFactor(
    double valueA,
    double valueB,
    double isovalue
)
{
    const double denominator =
        valueB - valueA;

    if (std::abs(denominator) < VALUE_EPSILON)
    {
        return 0.5;
    }

    return std::clamp(
        (isovalue - valueA) / denominator,
        0.0,
        1.0
    );
}

glm::vec3 calculateNormal(
    const glm::vec3& a,
    const glm::vec3& b,
    const glm::vec3& c
)
{
    glm::vec3 normal =
        glm::cross(
            b - a,
            c - a
        );

    const float length =
        glm::length(normal);

    if (length < 1.0e-12f)
    {
        return glm::vec3(0.0f);
    }

    return normal / length;
}

}

MolecularOrbitalIsosurface::Surface
MolecularOrbitalIsosurface::generate(
    const CartesianGrid& grid,
    const std::vector<double>& psi,
    double isovalue
)
{
    if (grid.getNx() < 2 ||
        grid.getNy() < 2 ||
        grid.getNz() < 2)
    {
        throw std::invalid_argument(
            "Molecular orbital grid must contain at least two points per dimension."
        );
    }

    if (psi.size() != grid.getSize())
    {
        throw std::invalid_argument(
            "Molecular orbital field size does not match the Cartesian grid."
        );
    }

    if (!std::isfinite(isovalue) ||
        isovalue <= 0.0)
    {
        throw std::invalid_argument(
            "Molecular orbital isovalue must be finite and greater than zero."
        );
    }

    Surface surface;

    /*
     * Cube vertices:
     *
     *             7 -------- 6
     *            /|         /|
     *           / |        / |
     *          4 -------- 5  |
     *          |  |        |  |
     *          |  3 -------- 2
     *          | /         | /
     *          |/          |/
     *          0 -------- 1
     *
     * The cube is divided using the diagonal
     *
     *              0 -------- 6
     *
     * producing six non-overlapping tetrahedra.
     */

    static constexpr int TETRAHEDRA[6][4] =
    {
        { 0, 1, 2, 6 },
        { 0, 2, 3, 6 },
        { 0, 3, 7, 6 },
        { 0, 7, 4, 6 },
        { 0, 4, 5, 6 },
        { 0, 5, 1, 6 }
    };

    for (std::size_t k = 0;
         k + 1 < grid.getNz();
         ++k)
    {
        for (std::size_t j = 0;
             j + 1 < grid.getNy();
             ++j)
        {
            for (std::size_t i = 0;
                 i + 1 < grid.getNx();
                 ++i)
            {
                GridPoint cube[8] =
                {
                    makeGridPoint(
                        grid,
                        psi,
                        i,
                        j,
                        k
                    ),

                    makeGridPoint(
                        grid,
                        psi,
                        i + 1,
                        j,
                        k
                    ),

                    makeGridPoint(
                        grid,
                        psi,
                        i + 1,
                        j + 1,
                        k
                    ),

                    makeGridPoint(
                        grid,
                        psi,
                        i,
                        j + 1,
                        k
                    ),

                    makeGridPoint(
                        grid,
                        psi,
                        i,
                        j,
                        k + 1
                    ),

                    makeGridPoint(
                        grid,
                        psi,
                        i + 1,
                        j,
                        k + 1
                    ),

                    makeGridPoint(
                        grid,
                        psi,
                        i + 1,
                        j + 1,
                        k + 1
                    ),

                    makeGridPoint(
                        grid,
                        psi,
                        i,
                        j + 1,
                        k + 1
                    )
                };

                for (const auto& tetrahedron : TETRAHEDRA)
                {
                    GridPoint tetrahedronPoints[4] =
                    {
                        cube[tetrahedron[0]],
                        cube[tetrahedron[1]],
                        cube[tetrahedron[2]],
                        cube[tetrahedron[3]]
                    };

                    polygonizeTetrahedron(
                        surface,
                        tetrahedronPoints,
                        isovalue
                    );

                    polygonizeTetrahedron(
                        surface,
                        tetrahedronPoints,
                        -isovalue
                    );
                }
            }
        }
    }

    return surface;
}

MolecularOrbitalIsosurface::GridPoint
MolecularOrbitalIsosurface::makeGridPoint(
    const CartesianGrid& grid,
    const std::vector<double>& psi,
    std::size_t i,
    std::size_t j,
    std::size_t k
)
{
    const std::size_t index =
        grid.getIndex(
            i,
            j,
            k
        );

    GridPoint point;

    point.position =
        glm::vec3(
            static_cast<float>(grid.getX(i)),
            static_cast<float>(grid.getY(j)),
            static_cast<float>(grid.getZ(k))
        );

    point.value =
        psi[index];

    return point;
}

glm::vec3
MolecularOrbitalIsosurface::interpolatePosition(
    const GridPoint& a,
    const GridPoint& b,
    double isovalue
)
{
    const double factor =
        interpolationFactor(
            a.value,
            b.value,
            isovalue
        );

    return a.position +
        static_cast<float>(factor) *
        (b.position - a.position);
}

void MolecularOrbitalIsosurface::addTriangle(
    Surface& surface,
    const glm::vec3& a,
    const glm::vec3& b,
    const glm::vec3& c,
    float phase
)
{
    const glm::vec3 normal =
        calculateNormal(
            a,
            b,
            c
        );

    const unsigned int firstIndex =
        static_cast<unsigned int>(
            surface.vertices.size()
        );

    surface.vertices.push_back(
        {
            a,
            normal,
            phase
        }
    );

    surface.vertices.push_back(
        {
            b,
            normal,
            phase
        }
    );

    surface.vertices.push_back(
        {
            c,
            normal,
            phase
        }
    );

    surface.indices.push_back(
        firstIndex
    );

    surface.indices.push_back(
        firstIndex + 1
    );

    surface.indices.push_back(
        firstIndex + 2
    );
}

void MolecularOrbitalIsosurface::polygonizeTetrahedron(
    Surface& surface,
    const GridPoint* points,
    double isovalue
)
{
    bool inside[4] =
    {
        false,
        false,
        false,
        false
    };

    int insideCount = 0;

    for (int i = 0; i < 4; ++i)
    {
        inside[i] =
            isInside(
                points[i].value,
                isovalue
            );

        if (inside[i])
        {
            ++insideCount;
        }
    }

    if (insideCount == 0 ||
        insideCount == 4)
    {
        return;
    }

    const float phase =
        isovalue > 0.0
            ? 1.0f
            : -1.0f;

    glm::vec3 intersections[6];

    int intersectionCount = 0;

    for (const TetrahedronEdge& edge :
         TETRAHEDRON_EDGES)
    {
        if (inside[edge.a] ==
            inside[edge.b])
        {
            continue;
        }

        intersections[intersectionCount++] =
            interpolatePosition(
                points[edge.a],
                points[edge.b],
                isovalue
            );
    }

    if (intersectionCount == 3)
    {
        addTriangle(
            surface,
            intersections[0],
            intersections[1],
            intersections[2],
            phase
        );

        return;
    }

    if (intersectionCount == 4)
    {
        addTriangle(
            surface,
            intersections[0],
            intersections[1],
            intersections[2],
            phase
        );

        addTriangle(
            surface,
            intersections[0],
            intersections[2],
            intersections[3],
            phase
        );
    }
}