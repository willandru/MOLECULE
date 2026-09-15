#include "DFTGrid.h"

#include <stdexcept>

DFTGrid::DFTGrid()
    : nx(0),
      ny(0),
      nz(0),
      spacing(0.0),
      origin(0.0),
      pointCount(0)
{
}

DFTGrid::DFTGrid(
    int nxValue,
    int nyValue,
    int nzValue,
    double spacingValue,
    const glm::dvec3& originValue
)
    : nx(0),
      ny(0),
      nz(0),
      spacing(0.0),
      origin(0.0),
      pointCount(0)
{
    initialize(
        nxValue,
        nyValue,
        nzValue,
        spacingValue,
        originValue
    );
}

void DFTGrid::initialize(
    int nxValue,
    int nyValue,
    int nzValue,
    double spacingValue,
    const glm::dvec3& originValue
)
{
    if (nxValue < 2 ||
        nyValue < 2 ||
        nzValue < 2)
    {
        throw std::invalid_argument(
            "DFTGrid requires at least 2 points per dimension."
        );
    }

    if (spacingValue <= 0.0)
    {
        throw std::invalid_argument(
            "DFTGrid spacing must be greater than zero."
        );
    }

    nx = nxValue;
    ny = nyValue;
    nz = nzValue;

    spacing = spacingValue;

    origin = originValue;

    pointCount =
        static_cast<std::size_t>(nx) *
        static_cast<std::size_t>(ny) *
        static_cast<std::size_t>(nz);
}

int DFTGrid::getNx() const
{
    return nx;
}

int DFTGrid::getNy() const
{
    return ny;
}

int DFTGrid::getNz() const
{
    return nz;
}

double DFTGrid::getSpacing() const
{
    return spacing;
}

const glm::dvec3& DFTGrid::getOrigin() const
{
    return origin;
}

std::size_t DFTGrid::getPointCount() const
{
    return pointCount;
}

double DFTGrid::getVolumeElement() const
{
    return spacing * spacing * spacing;
}

double DFTGrid::getVolume() const
{
    if (pointCount == 0)
        return 0.0;

    /*
        For the rectangular quadrature used by the DFT modules,
        each grid point represents a volume element dV.

        Therefore:

            V = Npoints * dV
    */

    return
        static_cast<double>(pointCount) *
        getVolumeElement();
}

std::size_t DFTGrid::getIndex(
    int x,
    int y,
    int z
) const
{
    if (!isInside(x, y, z))
    {
        throw std::out_of_range(
            "DFTGrid coordinates are outside the grid."
        );
    }

    /*
        Row-major layout:

            index = x + nx * (y + ny * z)
    */

    return
        static_cast<std::size_t>(x)
        +
        static_cast<std::size_t>(nx) *
        (
            static_cast<std::size_t>(y)
            +
            static_cast<std::size_t>(ny) *
            static_cast<std::size_t>(z)
        );
}

void DFTGrid::getCoordinates(
    std::size_t index,
    int& x,
    int& y,
    int& z
) const
{
    if (index >= pointCount)
    {
        throw std::out_of_range(
            "DFTGrid index is outside the grid."
        );
    }

    const std::size_t xy =
        static_cast<std::size_t>(nx) *
        static_cast<std::size_t>(ny);

    z =
        static_cast<int>(index / xy);

    const std::size_t remainder =
        index % xy;

    y =
        static_cast<int>(
            remainder /
            static_cast<std::size_t>(nx)
        );

    x =
        static_cast<int>(
            remainder %
            static_cast<std::size_t>(nx)
        );
}

glm::dvec3 DFTGrid::getPosition(
    int x,
    int y,
    int z
) const
{
    if (!isInside(x, y, z))
    {
        throw std::out_of_range(
            "DFTGrid coordinates are outside the grid."
        );
    }

    return origin +
           glm::dvec3(
               static_cast<double>(x) * spacing,
               static_cast<double>(y) * spacing,
               static_cast<double>(z) * spacing
           );
}

glm::dvec3 DFTGrid::getPosition(
    std::size_t index
) const
{
    int x;
    int y;
    int z;

    getCoordinates(
        index,
        x,
        y,
        z
    );

    return getPosition(x, y, z);
}

bool DFTGrid::isInside(
    int x,
    int y,
    int z
) const
{
    return
        x >= 0 && x < nx &&
        y >= 0 && y < ny &&
        z >= 0 && z < nz;
}