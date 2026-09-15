#include "PotentialEnergy.h"

#include <cmath>
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

PotentialEnergy::PotentialEnergy()
    :
    grid(nullptr)
{
}


PotentialEnergy::PotentialEnergy(
    const DFTGrid& grid
)
    :
    grid(nullptr)
{
    initialize(grid);
}


// ================================================================
// INICIALIZACIÓN
// ================================================================

void PotentialEnergy::initialize(
    const DFTGrid& grid
)
{
    if (grid.getPointCount() == 0)
    {
        throw std::invalid_argument(
            "PotentialEnergy: grid has no points."
        );
    }

    if (grid.getNx() < 3 ||
        grid.getNy() < 3 ||
        grid.getNz() < 3)
    {
        throw std::invalid_argument(
            "PotentialEnergy: grid dimensions must be at least 3 in every direction."
        );
    }

    if (!std::isfinite(grid.getSpacing()) ||
        grid.getSpacing() <= 0.0)
    {
        throw std::invalid_argument(
            "PotentialEnergy: grid spacing must be positive and finite."
        );
    }

    this->grid =
        &grid;
}


// ================================================================
// ENERGÍA ELECTRÓN - NÚCLEO
// ================================================================

double PotentialEnergy::calculateElectronNuclearEnergy(
    const DFTDensity& density,
    const std::vector<double>& externalPotential
) const
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "PotentialEnergy: not initialized."
        );
    }

    if (density.getGrid().getPointCount() !=
        grid->getPointCount())
    {
        throw std::invalid_argument(
            "PotentialEnergy: incompatible density grid."
        );
    }

    if (externalPotential.size() !=
        grid->getPointCount())
    {
        throw std::invalid_argument(
            "PotentialEnergy: incompatible external potential."
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
// ENERGÍA DE HARTREE
// ================================================================

double PotentialEnergy::calculateHartreeEnergy(
    const DFTDensity& density,
    const std::vector<double>& hartreePotential
) const
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "PotentialEnergy: not initialized."
        );
    }

    if (density.getGrid().getPointCount() !=
        grid->getPointCount())
    {
        throw std::invalid_argument(
            "PotentialEnergy: incompatible density grid."
        );
    }

    if (hartreePotential.size() !=
        grid->getPointCount())
    {
        throw std::invalid_argument(
            "PotentialEnergy: incompatible Hartree potential."
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
            hartreePotential[i] *
            dV;
    }

    return 0.5 * energy;
}


// ================================================================
// ENERGÍA DE INTERCAMBIO LDA
// ================================================================

double PotentialEnergy::calculateExchangeEnergy(
    const DFTDensity& density
) const
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "PotentialEnergy: not initialized."
        );
    }

    if (density.getGrid().getPointCount() !=
        grid->getPointCount())
    {
        throw std::invalid_argument(
            "PotentialEnergy: incompatible density grid."
        );
    }

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
            density.get(i);

        if (rho > MIN_DENSITY)
        {
            energy +=
                coefficient *
                std::pow(
                    rho,
                    4.0 / 3.0
                ) *
                dV;
        }
    }

    return energy;
}


// ================================================================
// REPULSIÓN NÚCLEO - NÚCLEO
// ================================================================

double PotentialEnergy::calculateNuclearRepulsionEnergy(
    const std::vector<int>& nuclearCharges,
    const std::vector<glm::dvec3>& nuclearPositions
) const
{
    if (nuclearCharges.size() !=
        nuclearPositions.size())
    {
        throw std::invalid_argument(
            "PotentialEnergy: nuclear charges and positions have different sizes."
        );
    }

    for (const int Z : nuclearCharges)
    {
        if (Z <= 0)
        {
            throw std::invalid_argument(
                "PotentialEnergy: nuclear charge must be positive."
            );
        }
    }

    double energy =
        0.0;

    for (std::size_t i = 0;
         i < nuclearCharges.size();
         ++i)
    {
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
                    "PotentialEnergy: coincident nuclear positions."
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