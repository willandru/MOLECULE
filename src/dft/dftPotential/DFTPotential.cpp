#include "DFTPotential.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace
{
    constexpr double PI =
        3.1415926535897932384626433832795;

    constexpr double DEFAULT_SOFTENING =
        0.15;

    constexpr double MIN_DENSITY =
        1.0e-14;

    constexpr double MIN_DISTANCE =
        1.0e-12;
}


// ================================================================
// CONSTRUCTORES
// ================================================================

DFTPotential::DFTPotential()
    :
    grid(nullptr),
    softening(DEFAULT_SOFTENING),
    multigridSolver()
{
}


DFTPotential::DFTPotential(
    const DFTGrid& grid
)
    :
    grid(nullptr),
    softening(DEFAULT_SOFTENING),
    multigridSolver()
{
    initialize(grid);
}


// ================================================================
// INICIALIZACIÓN
// ================================================================

void DFTPotential::initialize(
    const DFTGrid& grid
)
{
    if (grid.getPointCount() == 0)
    {
        throw std::invalid_argument(
            "DFTPotential: grid has no points."
        );
    }

    if (grid.getNx() < 3 ||
        grid.getNy() < 3 ||
        grid.getNz() < 3)
    {
        throw std::invalid_argument(
            "DFTPotential: grid dimensions must be at least 3 in every direction."
        );
    }

    if (!std::isfinite(grid.getSpacing()) ||
        grid.getSpacing() <= 0.0)
    {
        throw std::invalid_argument(
            "DFTPotential: grid spacing must be positive and finite."
        );
    }

    this->grid =
        &grid;

    const std::size_t count =
        grid.getPointCount();

    externalPotential.assign(
        count,
        0.0
    );

    hartreePotential.assign(
        count,
        0.0
    );

    exchangePotential.assign(
        count,
        0.0
    );

    multigridSolver.initialize(
        grid
    );
}


// ================================================================
// POTENCIAL EXTERNO
// ================================================================

void DFTPotential::calculateExternalPotential(
    const std::vector<int>& nuclearCharges,
    const std::vector<glm::dvec3>& nuclearPositions
)
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "DFTPotential: not initialized."
        );
    }

    if (nuclearCharges.size() !=
        nuclearPositions.size())
    {
        throw std::invalid_argument(
            "DFTPotential: nuclear charges and positions have different sizes."
        );
    }

    std::fill(
        externalPotential.begin(),
        externalPotential.end(),
        0.0
    );

    const std::size_t count =
        grid->getPointCount();

    for (std::size_t i = 0;
         i < count;
         ++i)
    {
        const glm::dvec3 r =
            grid->getPosition(i);

        double value =
            0.0;

        for (std::size_t nucleus = 0;
             nucleus < nuclearCharges.size();
             ++nucleus)
        {
            const int Z =
                nuclearCharges[nucleus];

            if (Z <= 0)
            {
                throw std::invalid_argument(
                    "DFTPotential: nuclear charge must be positive."
                );
            }

            const glm::dvec3 difference =
                r -
                nuclearPositions[nucleus];

            const double r2 =
                glm::dot(
                    difference,
                    difference
                );

            const double denominator =
                std::sqrt(
                    r2 +
                    softening * softening
                );

            value -=
                static_cast<double>(Z) /
                denominator;
        }

        externalPotential[i] =
            value;
    }
}


// ================================================================
// HARTREE
// ================================================================

void DFTPotential::calculateHartreePotential(
    const DFTDensity& density,
    int iterations
)
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "DFTPotential: not initialized."
        );
    }

    if (iterations <= 0)
    {
        throw std::invalid_argument(
            "DFTPotential: Hartree maximum cycles must be positive."
        );
    }

    if (density.getGrid().getPointCount() !=
        grid->getPointCount())
    {
        throw std::invalid_argument(
            "DFTPotential: incompatible density grid."
        );
    }

    /*
        Las condiciones de frontera pertenecen al
        problema original de Poisson.

        El MultigridSolver solamente modifica los
        puntos interiores.
    */

    applyBoundaryConditions(
        hartreePotential,
        density
    );

    multigridSolver.solve(
        density,
        hartreePotential,
        iterations
    );
}


// ================================================================
// LAPLACIANO
// ================================================================

double DFTPotential::calculateLaplacian(
    const std::vector<double>& potential,
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

    if (x <= 0 ||
        x >= nx - 1 ||
        y <= 0 ||
        y >= ny - 1 ||
        z <= 0 ||
        z >= nz - 1)
    {
        return 0.0;
    }

    const double h =
        grid->getSpacing();

    const double h2 =
        h * h;

    const std::size_t center =
        grid->getIndex(
            x,
            y,
            z
        );

    const double neighbors =
        potential[
            grid->getIndex(
                x + 1,
                y,
                z
            )
        ]
        +
        potential[
            grid->getIndex(
                x - 1,
                y,
                z
            )
        ]
        +
        potential[
            grid->getIndex(
                x,
                y + 1,
                z
            )
        ]
        +
        potential[
            grid->getIndex(
                x,
                y - 1,
                z
            )
        ]
        +
        potential[
            grid->getIndex(
                x,
                y,
                z + 1
            )
        ]
        +
        potential[
            grid->getIndex(
                x,
                y,
                z - 1
            )
        ];

    return
        (
            neighbors -
            6.0 * potential[center]
        )
        /
        h2;
}


// ================================================================
// CARGA ELECTRÓNICA TOTAL
// ================================================================

double DFTPotential::calculateDensityCharge(
    const DFTDensity& density
) const
{
    return density.calculateElectronCount();
}


// ================================================================
// CONDICIONES DE FRONTERA
// ================================================================

void DFTPotential::applyBoundaryConditions(
    std::vector<double>& potential,
    const DFTDensity& density
) const
{
    const int nx =
        grid->getNx();

    const int ny =
        grid->getNy();

    const int nz =
        grid->getNz();

    const double electronCount =
        calculateDensityCharge(
            density
        );

    if (!std::isfinite(electronCount) ||
        electronCount <= 0.0)
    {
        return;
    }

    const double dV =
        grid->getVolumeElement();

    glm::dvec3 center(
        0.0
    );

    double totalCharge =
        0.0;

    for (std::size_t i = 0;
         i < grid->getPointCount();
         ++i)
    {
        const double rho =
            density.get(i);

        if (!std::isfinite(rho))
        {
            throw std::runtime_error(
                "DFTPotential: density contains a non-finite value."
            );
        }

        center +=
            rho *
            grid->getPosition(i) *
            dV;

        totalCharge +=
            rho *
            dV;
    }

    if (totalCharge > MIN_DENSITY)
    {
        center /=
            totalCharge;
    }

    auto setBoundary =
        [&](int x, int y, int z)
        {
            const std::size_t index =
                grid->getIndex(
                    x,
                    y,
                    z
                );

            const glm::dvec3 difference =
                grid->getPosition(index) -
                center;

            const double distance =
                glm::length(
                    difference
                );

            if (distance > MIN_DISTANCE)
            {
                potential[index] =
                    electronCount /
                    distance;
            }
        };

    for (int x = 0;
         x < nx;
         ++x)
    {
        for (int y = 0;
             y < ny;
             ++y)
        {
            setBoundary(
                x,
                y,
                0
            );

            setBoundary(
                x,
                y,
                nz - 1
            );
        }
    }

    for (int x = 0;
         x < nx;
         ++x)
    {
        for (int z = 0;
             z < nz;
             ++z)
        {
            setBoundary(
                x,
                0,
                z
            );

            setBoundary(
                x,
                ny - 1,
                z
            );
        }
    }

    for (int y = 0;
         y < ny;
         ++y)
    {
        for (int z = 0;
             z < nz;
             ++z)
        {
            setBoundary(
                0,
                y,
                z
            );

            setBoundary(
                nx - 1,
                y,
                z
            );
        }
    }
}


// ================================================================
// INTERCAMBIO LDA
// ================================================================

void DFTPotential::calculateExchangePotential(
    const DFTDensity& density
)
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "DFTPotential: not initialized."
        );
    }

    if (density.getGrid().getPointCount() !=
        grid->getPointCount())
    {
        throw std::invalid_argument(
            "DFTPotential: incompatible density grid."
        );
    }

    const double coefficient =
        -std::cbrt(
            3.0 / PI
        );

    const std::size_t count =
        grid->getPointCount();

    for (std::size_t i = 0;
         i < count;
         ++i)
    {
        const double rho =
            std::max(
                density.get(i),
                0.0
            );

        if (rho <= MIN_DENSITY)
        {
            exchangePotential[i] =
                0.0;
        }
        else
        {
            exchangePotential[i] =
                coefficient *
                std::cbrt(
                    rho
                );
        }
    }
}


// ================================================================
// ACCESO
// ================================================================

double DFTPotential::getExternal(
    std::size_t index
) const
{
    return externalPotential.at(
        index
    );
}


double DFTPotential::getHartree(
    std::size_t index
) const
{
    return hartreePotential.at(
        index
    );
}


double DFTPotential::getExchange(
    std::size_t index
) const
{
    return exchangePotential.at(
        index
    );
}


double DFTPotential::getTotalElectronic(
    std::size_t index
) const
{
    return
        externalPotential.at(index)
        +
        hartreePotential.at(index)
        +
        exchangePotential.at(index);
}


const std::vector<double>&
DFTPotential::getExternalPotential() const
{
    return externalPotential;
}


const std::vector<double>&
DFTPotential::getHartreePotential() const
{
    return hartreePotential;
}


const std::vector<double>&
DFTPotential::getExchangePotential() const
{
    return exchangePotential;
}


// ================================================================
// DIAGNÓSTICO POISSON
// ================================================================

double DFTPotential::getPoissonResidual() const
{
    return multigridSolver.getResidual();
}


int DFTPotential::getPoissonCycles() const
{
    return multigridSolver.getCycles();
}


bool DFTPotential::hasPoissonConverged() const
{
    return multigridSolver.hasConverged();
}


// ================================================================
// ENERGÍA ELECTRÓN - NÚCLEO
// ================================================================

double DFTPotential::calculateElectronNuclearEnergy(
    const DFTDensity& density
) const
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "DFTPotential: not initialized."
        );
    }

    const double dV =
        grid->getVolumeElement();

    double energy =
        0.0;

    for (std::size_t i = 0;
         i < grid->getPointCount();
         ++i)
    {
        energy +=
            density.get(i) *
            externalPotential[i] *
            dV;
    }

    return energy;
}


// ================================================================
// ENERGÍA HARTREE
// ================================================================

double DFTPotential::calculateHartreeEnergy(
    const DFTDensity& density
) const
{
    const double dV =
        grid->getVolumeElement();

    double energy =
        0.0;

    for (std::size_t i = 0;
         i < grid->getPointCount();
         ++i)
    {
        energy +=
            density.get(i) *
            hartreePotential[i] *
            dV;
    }

    return
        0.5 *
        energy;
}


// ================================================================
// ENERGÍA DE INTERCAMBIO
// ================================================================

double DFTPotential::calculateExchangeEnergy(
    const DFTDensity& density
) const
{
    const double coefficient =
        -0.75 *
        std::cbrt(
            3.0 / PI
        );

    const double dV =
        grid->getVolumeElement();

    double energy =
        0.0;

    for (std::size_t i = 0;
         i < grid->getPointCount();
         ++i)
    {
        const double rho =
            std::max(
                density.get(i),
                0.0
            );

        if (rho > MIN_DENSITY)
        {
            energy +=
                coefficient *
                std::pow(
                    rho,
                    4.0 / 3.0
                )
                *
                dV;
        }
    }

    return energy;
}


// ================================================================
// REPULSIÓN NÚCLEO - NÚCLEO
// ================================================================

double DFTPotential::calculateNuclearRepulsionEnergy(
    const std::vector<int>& nuclearCharges,
    const std::vector<glm::dvec3>& nuclearPositions
) const
{
    if (nuclearCharges.size() !=
        nuclearPositions.size())
    {
        throw std::invalid_argument(
            "DFTPotential: nuclear charges and positions have different sizes."
        );
    }

    double energy =
        0.0;

    for (std::size_t i = 0;
         i < nuclearCharges.size();
         ++i)
    {
        if (nuclearCharges[i] <= 0)
        {
            throw std::invalid_argument(
                "DFTPotential: nuclear charge must be positive."
            );
        }

        for (std::size_t j = i + 1;
             j < nuclearCharges.size();
             ++j)
        {
            const glm::dvec3 difference =
                nuclearPositions[i] -
                nuclearPositions[j];

            const double distance =
                glm::length(
                    difference
                );

            if (distance <= MIN_DISTANCE)
            {
                throw std::runtime_error(
                    "DFTPotential: nuclei overlap."
                );
            }

            energy +=
                static_cast<double>(
                    nuclearCharges[i]
                )
                *
                static_cast<double>(
                    nuclearCharges[j]
                )
                /
                distance;
        }
    }

    return energy;
}


// ================================================================
// FUERZAS NUCLEARES
// ================================================================

std::vector<glm::dvec3>
DFTPotential::calculateNuclearForces(
    const DFTDensity& density,
    const std::vector<int>& nuclearCharges,
    const std::vector<glm::dvec3>& nuclearPositions
) const
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "DFTPotential: not initialized."
        );
    }

    if (nuclearCharges.size() !=
        nuclearPositions.size())
    {
        throw std::invalid_argument(
            "DFTPotential: nuclear charges and positions have different sizes."
        );
    }

    std::vector<glm::dvec3> forces(
        nuclearPositions.size(),
        glm::dvec3(0.0)
    );

    const double dV =
        grid->getVolumeElement();

    // ------------------------------------------------------------
    // Fuerza electrón-núcleo
    //
    // Consistente con el potencial suavizado
    // utilizado actualmente.
    // ------------------------------------------------------------

    for (std::size_t nucleus = 0;
         nucleus < nuclearPositions.size();
         ++nucleus)
    {
        const double Z =
            static_cast<double>(
                nuclearCharges[nucleus]
            );

        glm::dvec3 electronicForce(
            0.0
        );

        for (std::size_t i = 0;
             i < grid->getPointCount();
             ++i)
        {
            const glm::dvec3 difference =
                grid->getPosition(i) -
                nuclearPositions[nucleus];

            const double r2 =
                glm::dot(
                    difference,
                    difference
                );

            const double denominator =
                std::pow(
                    r2 +
                    softening * softening,
                    1.5
                );

            electronicForce +=
                Z *
                density.get(i) *
                difference /
                denominator *
                dV;
        }

        forces[nucleus] -=
            electronicForce;
    }

    // ------------------------------------------------------------
    // Fuerza núcleo-núcleo
    // ------------------------------------------------------------

    for (std::size_t i = 0;
         i < nuclearPositions.size();
         ++i)
    {
        for (std::size_t j = i + 1;
             j < nuclearPositions.size();
             ++j)
        {
            const glm::dvec3 difference =
                nuclearPositions[i] -
                nuclearPositions[j];

            const double r2 =
                glm::dot(
                    difference,
                    difference
                );

            const double distance =
                std::sqrt(
                    r2
                );

            if (distance <= MIN_DISTANCE)
            {
                throw std::runtime_error(
                    "DFTPotential: nuclei overlap."
                );
            }

            const glm::dvec3 force =
                static_cast<double>(
                    nuclearCharges[i]
                )
                *
                static_cast<double>(
                    nuclearCharges[j]
                )
                *
                difference
                /
                (
                    distance *
                    distance *
                    distance
                );

            forces[i] +=
                force;

            forces[j] -=
                force;
        }
    }

    return forces;
}


// ================================================================
// SOFTENING
// ================================================================

void DFTPotential::setSoftening(
    double value
)
{
    if (!std::isfinite(value) ||
        value <= 0.0)
    {
        throw std::invalid_argument(
            "DFTPotential: softening must be positive and finite."
        );
    }

    softening =
        value;
}


double DFTPotential::getSoftening() const
{
    return softening;
}


// ================================================================
// GRID
// ================================================================

const DFTGrid&
DFTPotential::getGrid() const
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "DFTPotential: not initialized."
        );
    }

    return *grid;
}