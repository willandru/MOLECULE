#include "ExternalPotential.h"

#include <cmath>
#include <stdexcept>

namespace
{
    constexpr double DEFAULT_SOFTENING =
        0.15;
}


// ================================================================
// CONSTRUCTORES
// ================================================================

ExternalPotential::ExternalPotential()
    :
    grid(nullptr),
    potential(),
    softening(DEFAULT_SOFTENING)
{
}


ExternalPotential::ExternalPotential(
    const DFTGrid& grid
)
    :
    grid(nullptr),
    potential(),
    softening(DEFAULT_SOFTENING)
{
    initialize(grid);
}


// ================================================================
// INICIALIZACIÓN
// ================================================================

void ExternalPotential::initialize(
    const DFTGrid& grid
)
{
    if (grid.getPointCount() == 0)
    {
        throw std::invalid_argument(
            "ExternalPotential: grid has no points."
        );
    }

    if (grid.getNx() < 3 ||
        grid.getNy() < 3 ||
        grid.getNz() < 3)
    {
        throw std::invalid_argument(
            "ExternalPotential: grid dimensions must be at least 3 in every direction."
        );
    }

    if (!std::isfinite(grid.getSpacing()) ||
        grid.getSpacing() <= 0.0)
    {
        throw std::invalid_argument(
            "ExternalPotential: grid spacing must be positive and finite."
        );
    }

    this->grid =
        &grid;

    potential.assign(
        grid.getPointCount(),
        0.0
    );
}


// ================================================================
// CÁLCULO
// ================================================================

void ExternalPotential::calculate(
    const std::vector<int>& nuclearCharges,
    const std::vector<glm::dvec3>& nuclearPositions
)
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "ExternalPotential: not initialized."
        );
    }

    if (nuclearCharges.size() !=
        nuclearPositions.size())
    {
        throw std::invalid_argument(
            "ExternalPotential: nuclear charges and positions have different sizes."
        );
    }

    for (const int Z : nuclearCharges)
    {
        if (Z <= 0)
        {
            throw std::invalid_argument(
                "ExternalPotential: nuclear charge must be positive."
            );
        }
    }

    std::fill(
        potential.begin(),
        potential.end(),
        0.0
    );

    for (std::size_t i = 0;
         i < grid->getPointCount();
         ++i)
    {
        const glm::dvec3 position =
            grid->getPosition(i);

        for (std::size_t nucleus = 0;
             nucleus < nuclearCharges.size();
             ++nucleus)
        {
            const glm::dvec3 difference =
                position -
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

            potential[i] -=
                static_cast<double>(
                    nuclearCharges[nucleus]
                )
                /
                denominator;
        }
    }
}


// ================================================================
// ACCESO
// ================================================================

double ExternalPotential::get(
    std::size_t index
) const
{
    return potential.at(
        index
    );
}


const std::vector<double>&
ExternalPotential::getPotential() const
{
    return potential;
}


// ================================================================
// SOFTENING
// ================================================================

void ExternalPotential::setSoftening(
    double value
)
{
    if (!std::isfinite(value) ||
        value <= 0.0)
    {
        throw std::invalid_argument(
            "ExternalPotential: softening must be positive and finite."
        );
    }

    softening =
        value;
}


double ExternalPotential::getSoftening() const
{
    return softening;
}