#pragma once

#include <vector>
#include <cstddef>

#include "DFTGrid.h"
#include "DFTDensity.h"

class ExchangePotential
{
public:

    ExchangePotential();

    explicit ExchangePotential(
        const DFTGrid& grid
    );

    void initialize(
        const DFTGrid& grid
    );

    void calculate(
        const DFTDensity& density
    );

    double get(
        std::size_t index
    ) const;

    const std::vector<double>& getPotential() const;

private:

    const DFTGrid* grid;

    std::vector<double> potential;
};