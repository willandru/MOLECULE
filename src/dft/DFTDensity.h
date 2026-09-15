#pragma once

#include <vector>
#include <cstddef>

#include "DFTGrid.h"

class DFTDensity
{
public:

    DFTDensity();

    explicit DFTDensity(
        const DFTGrid& grid
    );

    void initialize(
        const DFTGrid& grid
    );

    // Acceso a la densidad
    double get(std::size_t index) const;

    double& get(std::size_t index);

    void set(
        std::size_t index,
        double value
    );

    // Acceso directo al vector
    const std::vector<double>& getValues() const;

    std::vector<double>& getValues();

    // Tamaño
    std::size_t size() const;

    // Operaciones fundamentales
    void clear();

    void fill(double value);

    void normalize(
        double electronCount
    );

    double integrate() const;

    double calculateElectronCount() const;

    double getMaximum() const;

    double getMinimum() const;

    // Diferencia entre dos densidades
    double calculateDifferenceNorm(
        const DFTDensity& other
    ) const;

    double calculateMaximumDifference(
        const DFTDensity& other
    ) const;

    // Operaciones útiles para SCF
    void mix(
        const DFTDensity& other,
        double mixing
    );

    // Estadísticas
    double calculateMean() const;

    // Grid asociado
    const DFTGrid& getGrid() const;

private:

    const DFTGrid* grid;

    std::vector<double> density;
};