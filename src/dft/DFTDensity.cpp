#include "DFTDensity.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

DFTDensity::DFTDensity()
    : grid(nullptr)
{
}

DFTDensity::DFTDensity(
    const DFTGrid& gridReference
)
    : grid(nullptr)
{
    initialize(gridReference);
}

void DFTDensity::initialize(
    const DFTGrid& gridReference
)
{
    grid = &gridReference;

    density.assign(
        grid->getPointCount(),
        0.0
    );
}

double DFTDensity::get(
    std::size_t index
) const
{
    if (index >= density.size())
    {
        throw std::out_of_range(
            "DFTDensity index is outside the density array."
        );
    }

    return density[index];
}

double& DFTDensity::get(
    std::size_t index
)
{
    if (index >= density.size())
    {
        throw std::out_of_range(
            "DFTDensity index is outside the density array."
        );
    }

    return density[index];
}

void DFTDensity::set(
    std::size_t index,
    double value
)
{
    if (index >= density.size())
    {
        throw std::out_of_range(
            "DFTDensity index is outside the density array."
        );
    }

    if (!std::isfinite(value))
    {
        throw std::invalid_argument(
            "DFTDensity cannot contain NaN or infinity."
        );
    }

    /*
        Electron density cannot be negative.

        Small negative values can sometimes appear from numerical
        noise when a density is constructed from orbitals. We do
        not silently accept them.
    */
    if (value < 0.0)
    {
        throw std::invalid_argument(
            "DFTDensity cannot contain negative values."
        );
    }

    density[index] = value;
}

const std::vector<double>&
DFTDensity::getValues() const
{
    return density;
}

std::vector<double>&
DFTDensity::getValues()
{
    return density;
}

std::size_t DFTDensity::size() const
{
    return density.size();
}

void DFTDensity::clear()
{
    std::fill(
        density.begin(),
        density.end(),
        0.0
    );
}

void DFTDensity::fill(
    double value
)
{
    if (!std::isfinite(value))
    {
        throw std::invalid_argument(
            "DFTDensity cannot be filled with NaN or infinity."
        );
    }

    if (value < 0.0)
    {
        throw std::invalid_argument(
            "DFTDensity cannot be negative."
        );
    }

    std::fill(
        density.begin(),
        density.end(),
        value
    );
}

void DFTDensity::normalize(
    double electronCount
)
{
    if (electronCount < 0.0)
    {
        throw std::invalid_argument(
            "Electron count cannot be negative."
        );
    }

    if (density.empty())
    {
        throw std::runtime_error(
            "Cannot normalize an empty density."
        );
    }

    const double currentElectronCount =
        integrate();

    if (currentElectronCount <= 0.0)
    {
        throw std::runtime_error(
            "Cannot normalize a density with zero integral."
        );
    }

    const double scale =
        electronCount /
        currentElectronCount;

    for (double& value : density)
        value *= scale;
}

double DFTDensity::integrate() const
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "DFTDensity has no associated grid."
        );
    }

    double sum = 0.0;

    for (const double value : density)
        sum += value;

    return
        sum *
        grid->getVolumeElement();
}

double DFTDensity::calculateElectronCount() const
{
    return integrate();
}

double DFTDensity::getMaximum() const
{
    if (density.empty())
        return 0.0;

    return *std::max_element(
        density.begin(),
        density.end()
    );
}

double DFTDensity::getMinimum() const
{
    if (density.empty())
        return 0.0;

    return *std::min_element(
        density.begin(),
        density.end()
    );
}

double DFTDensity::calculateDifferenceNorm(
    const DFTDensity& other
) const
{
    if (grid == nullptr ||
        other.grid == nullptr)
    {
        throw std::runtime_error(
            "Both densities must have associated grids."
        );
    }

    if (density.size() != other.density.size())
    {
        throw std::invalid_argument(
            "Cannot compare densities with different sizes."
        );
    }

    /*
        L2 norm of the density difference:

            ||ρ₁-ρ₂||
            =
            sqrt(
                ∫ |ρ₁-ρ₂|² dr
            )
    */

    double sum = 0.0;

    for (std::size_t i = 0;
         i < density.size();
         ++i)
    {
        const double difference =
            density[i] -
            other.density[i];

        sum +=
            difference *
            difference;
    }

    return std::sqrt(
        sum *
        grid->getVolumeElement()
    );
}

double DFTDensity::calculateMaximumDifference(
    const DFTDensity& other
) const
{
    if (density.size() != other.density.size())
    {
        throw std::invalid_argument(
            "Cannot compare densities with different sizes."
        );
    }

    double maximumDifference = 0.0;

    for (std::size_t i = 0;
         i < density.size();
         ++i)
    {
        maximumDifference =
            std::max(
                maximumDifference,
                std::abs(
                    density[i] -
                    other.density[i]
                )
            );
    }

    return maximumDifference;
}

void DFTDensity::mix(
    const DFTDensity& other,
    double mixing
)
{
    if (density.size() != other.density.size())
    {
        throw std::invalid_argument(
            "Cannot mix densities with different sizes."
        );
    }

    if (mixing < 0.0 ||
        mixing > 1.0)
    {
        throw std::invalid_argument(
            "Density mixing parameter must be between 0 and 1."
        );
    }

    /*
        Linear mixing:

            ρ_new =
                (1-α)ρ_old +
                α ρ_calculated

        where α = mixing.
    */

    for (std::size_t i = 0;
         i < density.size();
         ++i)
    {
        density[i] =
            (1.0 - mixing) * density[i]
            +
            mixing * other.density[i];

        /*
            Protect against tiny negative values produced by
            floating-point arithmetic.
        */
        if (density[i] < 0.0 &&
            density[i] > -1.0e-14)
        {
            density[i] = 0.0;
        }

        if (density[i] < 0.0)
        {
            throw std::runtime_error(
                "Density mixing produced a negative density."
            );
        }
    }
}

double DFTDensity::calculateMean() const
{
    if (density.empty())
        return 0.0;

    double sum = 0.0;

    for (double value : density)
        sum += value;

    return
        sum /
        static_cast<double>(density.size());
}

const DFTGrid& DFTDensity::getGrid() const
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "DFTDensity has no associated grid."
        );
    }

    return *grid;
}