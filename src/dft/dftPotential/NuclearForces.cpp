#include "NuclearForces.h"

#include <cmath>
#include <stdexcept>


// ================================================================
// CONSTRUCTORES
// ================================================================

NuclearForces::NuclearForces()
    :
    grid(nullptr)
{
}


NuclearForces::NuclearForces(
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

void NuclearForces::initialize(
    const DFTGrid& grid
)
{
    if (grid.getPointCount() == 0)
    {
        throw std::invalid_argument(
            "NuclearForces: grid has no points."
        );
    }

    if (grid.getNx() < 3 ||
        grid.getNy() < 3 ||
        grid.getNz() < 3)
    {
        throw std::invalid_argument(
            "NuclearForces: grid dimensions must be at least 3 in every direction."
        );
    }

    if (!std::isfinite(grid.getSpacing()) ||
        grid.getSpacing() <= 0.0)
    {
        throw std::invalid_argument(
            "NuclearForces: grid spacing must be positive and finite."
        );
    }

    this->grid =
        &grid;
}


// ================================================================
// CÁLCULO DE FUERZAS NUCLEARES
// ================================================================

std::vector<glm::dvec3> NuclearForces::calculate(
    const DFTDensity& density,
    const std::vector<int>& nuclearCharges,
    const std::vector<glm::dvec3>& nuclearPositions,
    double softening
) const
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "NuclearForces: not initialized."
        );
    }

    if (density.getGrid().getPointCount() !=
        grid->getPointCount())
    {
        throw std::invalid_argument(
            "NuclearForces: incompatible density grid."
        );
    }

    if (nuclearCharges.size() !=
        nuclearPositions.size())
    {
        throw std::invalid_argument(
            "NuclearForces: nuclear charges and positions have different sizes."
        );
    }

    if (!std::isfinite(softening) ||
        softening <= 0.0)
    {
        throw std::invalid_argument(
            "NuclearForces: softening must be positive and finite."
        );
    }

    for (const int Z : nuclearCharges)
    {
        if (Z <= 0)
        {
            throw std::invalid_argument(
                "NuclearForces: nuclear charge must be positive."
            );
        }
    }

    const std::size_t nuclearCount =
        nuclearCharges.size();

    std::vector<glm::dvec3> forces(
        nuclearCount,
        glm::dvec3(0.0)
    );

    const double dV =
        grid->getVolumeElement();

    // ------------------------------------------------------------
    // Fuerza electrón - núcleo
    // ------------------------------------------------------------

    for (std::size_t nucleus = 0;
         nucleus < nuclearCount;
         ++nucleus)
    {
        glm::dvec3 force(
            0.0
        );

        const double Z =
            static_cast<double>(
                nuclearCharges[nucleus]
            );

        for (std::size_t i = 0;
             i < grid->getPointCount();
             ++i)
        {
            const double rho =
                density.get(i);

            if (!std::isfinite(rho))
            {
                throw std::runtime_error(
                    "NuclearForces: density contains a non-finite value."
                );
            }

            if (rho <= 0.0)
            {
                continue;
            }

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

            force +=
                Z *
                rho *
                difference /
                denominator *
                dV;
        }

        forces[nucleus] +=
            force;
    }

    // ------------------------------------------------------------
    // Fuerza núcleo - núcleo
    // ------------------------------------------------------------

    for (std::size_t i = 0;
         i < nuclearCount;
         ++i)
    {
        for (std::size_t j = i + 1;
             j < nuclearCount;
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

            if (distance <= 0.0)
            {
                throw std::runtime_error(
                    "NuclearForces: coincident nuclear positions."
                );
            }

            const double coefficient =
                static_cast<double>(
                    nuclearCharges[i]
                )
                *
                static_cast<double>(
                    nuclearCharges[j]
                )
                /
                (r2 * distance);

            const glm::dvec3 force =
                coefficient *
                difference;

            forces[i] +=
                force;

            forces[j] -=
                force;
        }
    }

    return forces;
}