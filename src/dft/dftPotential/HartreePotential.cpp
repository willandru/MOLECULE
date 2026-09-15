#include "HartreePotential.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>


namespace
{
    constexpr double PI =
        3.1415926535897932384626433832795;

    constexpr double MIN_DENSITY =
        1.0e-14;

    constexpr double MIN_DISTANCE =
        1.0e-12;
}


// ================================================================
// CONSTRUCTORES
// ================================================================

HartreePotential::HartreePotential()
    :
    grid(nullptr),
    potential(),
    multigridSolver()
{
}


HartreePotential::HartreePotential(
    const DFTGrid& grid
)
    :
    grid(nullptr),
    potential(),
    multigridSolver()
{
    initialize(grid);
}


// ================================================================
// INICIALIZACIÓN
// ================================================================

void HartreePotential::initialize(
    const DFTGrid& grid
)
{
    if (grid.getPointCount() == 0)
    {
        throw std::invalid_argument(
            "HartreePotential: grid has no points."
        );
    }

    if (grid.getNx() < 3 ||
        grid.getNy() < 3 ||
        grid.getNz() < 3)
    {
        throw std::invalid_argument(
            "HartreePotential: grid dimensions must be at least 3 in every direction."
        );
    }

    if (!std::isfinite(grid.getSpacing()) ||
        grid.getSpacing() <= 0.0)
    {
        throw std::invalid_argument(
            "HartreePotential: grid spacing must be positive and finite."
        );
    }

    this->grid =
        &grid;

    potential.assign(
        grid.getPointCount(),
        0.0
    );

    multigridSolver.initialize(
        grid
    );
}


// ================================================================
// CÁLCULO DEL POTENCIAL DE HARTREE
// ================================================================

void HartreePotential::calculate(
    const DFTDensity& density,
    int iterations
)
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "HartreePotential: not initialized."
        );
    }

    if (iterations <= 0)
    {
        throw std::invalid_argument(
            "HartreePotential: maximum cycles must be positive."
        );
    }

    if (density.getGrid().getPointCount() !=
        grid->getPointCount())
    {
        throw std::invalid_argument(
            "HartreePotential: incompatible density grid."
        );
    }

    applyBoundaryConditions(
        density
    );

    multigridSolver.solve(
        density,
        potential,
        iterations
    );

    std::cout
        << "[HartreePotential] solve"
        << " | residual = "
        << multigridSolver.getResidual()
        << " | cycles = "
        << multigridSolver.getCycles()
        << " | converged = "
        << (
            multigridSolver.hasConverged()
            ? "yes"
            : "no"
        )
        << '\n';
}


// ================================================================
// CONDICIONES DE FRONTERA MULTIPOLARES
// ================================================================
//
// Para un sistema electrónico localizado:
//
//     V_H(r) = integral rho(r') / |r-r'| dr'
//
// Lejos de la distribución:
//
//     V_H(r) = Q/r
//            + p.r/r^3
//            + términos de orden superior
//
// donde:
//
//     Q = integral rho(r') dr'
//
// es la carga electrónica total.
//
// Para un átomo, la distribución es aproximadamente esférica
// y el término monopolar domina.
//
// Para una molécula, los términos dipolar y cuadrupolar permiten
// representar una distribución electrónica no esférica.
//
// Se utilizan coordenadas relativas al centro de carga electrónica.
// ================================================================

void HartreePotential::applyBoundaryConditions(
    const DFTDensity& density
)
{
    const int nx =
        grid->getNx();

    const int ny =
        grid->getNy();

    const int nz =
        grid->getNz();

    const double dV =
        grid->getVolumeElement();

    const std::size_t pointCount =
        grid->getPointCount();

    // ------------------------------------------------------------
    // 1. Carga electrónica total
    // ------------------------------------------------------------

    double totalCharge =
        0.0;

    glm::dvec3 center(
        0.0
    );

    for (std::size_t i = 0;
         i < pointCount;
         ++i)
    {
        const double rho =
            density.get(i);

        if (!std::isfinite(rho))
        {
            throw std::runtime_error(
                "HartreePotential: density contains a non-finite value."
            );
        }

        if (rho < 0.0)
        {
            throw std::runtime_error(
                "HartreePotential: density contains a negative value."
            );
        }

        totalCharge +=
            rho *
            dV;

        center +=
            rho *
            grid->getPosition(i) *
            dV;
    }

    if (!std::isfinite(totalCharge) ||
        totalCharge <= MIN_DENSITY)
    {
        std::fill(
            potential.begin(),
            potential.end(),
            0.0
        );

        std::cout
            << "[HartreePotential] boundary"
            << " | electron count = 0"
            << " | potential boundary set to zero"
            << '\n';

        return;
    }

    center /=
        totalCharge;

    // ------------------------------------------------------------
    // 2. Momento dipolar
    // ------------------------------------------------------------

    glm::dvec3 dipole(
        0.0
    );

    // ------------------------------------------------------------
    // 3. Momento cuadrupolar
    //
    // Tensor cartesiano:
    //
    // Q_ij =
    // integral rho *
    // (3 r_i r_j - r^2 delta_ij) dr
    // ------------------------------------------------------------

    double Qxx =
        0.0;

    double Qxy =
        0.0;

    double Qxz =
        0.0;

    double Qyy =
        0.0;

    double Qyz =
        0.0;

    double Qzz =
        0.0;

    for (std::size_t i = 0;
         i < pointCount;
         ++i)
    {
        const double rho =
            density.get(i);

        if (rho <= 0.0)
        {
            continue;
        }

        const glm::dvec3 relative =
            grid->getPosition(i) -
            center;

        const double x =
            relative.x;

        const double y =
            relative.y;

        const double z =
            relative.z;

        const double r2 =
            x * x +
            y * y +
            z * z;

        dipole +=
            rho *
            relative *
            dV;

        Qxx +=
            rho *
            (
                3.0 * x * x -
                r2
            ) *
            dV;

        Qxy +=
            rho *
            (
                3.0 * x * y
            ) *
            dV;

        Qxz +=
            rho *
            (
                3.0 * x * z
            ) *
            dV;

        Qyy +=
            rho *
            (
                3.0 * y * y -
                r2
            ) *
            dV;

        Qyz +=
            rho *
            (
                3.0 * y * z
            ) *
            dV;

        Qzz +=
            rho *
            (
                3.0 * z * z -
                r2
            ) *
            dV;
    }

    // ------------------------------------------------------------
    // 4. Potencial multipolar en la frontera
    // ------------------------------------------------------------

    double minimumBoundaryPotential =
        std::numeric_limits<double>::max();

    double maximumBoundaryPotential =
        -std::numeric_limits<double>::max();

    std::size_t boundaryPointCount =
        0;

    for (int z = 0;
         z < nz;
         ++z)
    {
        for (int y = 0;
             y < ny;
             ++y)
        {
            for (int x = 0;
                 x < nx;
                 ++x)
            {
                if (!isBoundaryPoint(
                        x,
                        y,
                        z))
                {
                    continue;
                }

                const std::size_t index =
                    grid->getIndex(
                        x,
                        y,
                        z
                    );

                const glm::dvec3 relative =
                    grid->getPosition(index) -
                    center;

                const double rx =
                    relative.x;

                const double ry =
                    relative.y;

                const double rz =
                    relative.z;

                const double r2 =
                    rx * rx +
                    ry * ry +
                    rz * rz;

                const double r =
                    std::sqrt(
                        r2
                    );

                if (!std::isfinite(r) ||
                    r <= MIN_DISTANCE)
                {
                    throw std::runtime_error(
                        "HartreePotential: invalid boundary distance."
                    );
                }

                // ------------------------------------------------
                // Monopolo
                // ------------------------------------------------

                const double monopole =
                    totalCharge /
                    r;

                // ------------------------------------------------
                // Dipolo
                //
                // p.r / r^3
                // ------------------------------------------------

                const double dipoleTerm =
                    glm::dot(
                        dipole,
                        relative
                    )
                    /
                    (r2 * r);

                // ------------------------------------------------
                // Cuadrupolo
                //
                // 1/2 * r_i Q_ij r_j / r^5
                // ------------------------------------------------

                const double quadraticForm =
                    Qxx * rx * rx
                    +
                    2.0 * Qxy * rx * ry
                    +
                    2.0 * Qxz * rx * rz
                    +
                    Qyy * ry * ry
                    +
                    2.0 * Qyz * ry * rz
                    +
                    Qzz * rz * rz;

                const double quadrupoleTerm =
                    0.5 *
                    quadraticForm
                    /
                    (r2 * r2 * r);

                const double boundaryPotential =
                    monopole
                    +
                    dipoleTerm
                    +
                    quadrupoleTerm;

                if (!std::isfinite(
                        boundaryPotential))
                {
                    throw std::runtime_error(
                        "HartreePotential: non-finite boundary potential."
                    );
                }

                potential[index] =
                    boundaryPotential;

                minimumBoundaryPotential =
                    std::min(
                        minimumBoundaryPotential,
                        boundaryPotential
                    );

                maximumBoundaryPotential =
                    std::max(
                        maximumBoundaryPotential,
                        boundaryPotential
                    );

                ++boundaryPointCount;
            }
        }
    }

    std::cout
        << "[HartreePotential] boundary"
        << " | electron count = "
        << totalCharge
        << " | center = ("
        << center.x
        << ", "
        << center.y
        << ", "
        << center.z
        << ")"
        << " | boundary points = "
        << boundaryPointCount
        << " | min VH = "
        << minimumBoundaryPotential
        << " | max VH = "
        << maximumBoundaryPotential
        << '\n';
}


// ================================================================
// IDENTIFICACIÓN DE FRONTERA
// ================================================================

bool HartreePotential::isBoundaryPoint(
    int x,
    int y,
    int z
) const
{
    const int nx =
        grid->getNx();

    const int ny =
        grid->getNy();

    const int nz =
        grid->getNz();

    return
        x == 0 ||
        x == nx - 1 ||
        y == 0 ||
        y == ny - 1 ||
        z == 0 ||
        z == nz - 1;
}


// ================================================================
// ACCESO
// ================================================================

double HartreePotential::get(
    std::size_t index
) const
{
    return potential.at(
        index
    );
}


const std::vector<double>&
HartreePotential::getPotential() const
{
    return potential;
}


// ================================================================
// DIAGNÓSTICO DEL SOLVER
// ================================================================

double HartreePotential::getResidual() const
{
    return multigridSolver.getResidual();
}


int HartreePotential::getCycles() const
{
    return multigridSolver.getCycles();
}


bool HartreePotential::hasConverged() const
{
    return multigridSolver.hasConverged();
}