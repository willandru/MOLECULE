#include "ExchangePotential.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace
{
    constexpr double PI =
        3.1415926535897932384626433832795;

    constexpr double MIN_DENSITY =
        1.0e-14;
}


// ================================================================
// CONSTRUCTORES
// ================================================================

ExchangePotential::ExchangePotential()
    :
    grid(nullptr),
    potential()
{
}


ExchangePotential::ExchangePotential(
    const DFTGrid& grid
)
    :
    grid(nullptr),
    potential()
{
    initialize(grid);
}


// ================================================================
// INICIALIZACIÓN
// ================================================================

void ExchangePotential::initialize(
    const DFTGrid& grid
)
{
    if (grid.getPointCount() == 0)
    {
        throw std::invalid_argument(
            "ExchangePotential: grid has no points."
        );
    }

    if (grid.getNx() < 3 ||
        grid.getNy() < 3 ||
        grid.getNz() < 3)
    {
        throw std::invalid_argument(
            "ExchangePotential: grid dimensions must be at least 3 in every direction."
        );
    }

    if (!std::isfinite(grid.getSpacing()) ||
        grid.getSpacing() <= 0.0)
    {
        throw std::invalid_argument(
            "ExchangePotential: grid spacing must be positive and finite."
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
// CÁLCULO DEL POTENCIAL DE INTERCAMBIO LDA
// ================================================================

void ExchangePotential::calculate(
    const DFTDensity& density
)
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "ExchangePotential: not initialized."
        );
    }

    if (density.getGrid().getPointCount() !=
        grid->getPointCount())
    {
        throw std::invalid_argument(
            "ExchangePotential: incompatible density grid."
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
            potential[i] =
                0.0;
        }
        else
        {
            potential[i] =
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

double ExchangePotential::get(
    std::size_t index
) const
{
    return potential.at(
        index
    );
}


const std::vector<double>&
ExchangePotential::getPotential() const
{
    return potential;
}