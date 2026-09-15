#pragma once

#include <vector>
#include <cstddef>

#include "DFTGrid.h"
#include "DFTDensity.h"

class MultigridSolver
{
public:

    MultigridSolver();

    explicit MultigridSolver(
        const DFTGrid& grid
    );

    void initialize(
        const DFTGrid& grid
    );

    void solve(
        const DFTDensity& density,
        std::vector<double>& potential,
        int maxCycles
    );

    double getResidual() const;

    int getCycles() const;

    bool hasConverged() const;

private:

    const DFTGrid* grid;

    double residual;

    int cycles;

    bool converged;

    // ------------------------------------------------------------
    // V-cycle
    // ------------------------------------------------------------

    void vCycle(
        std::vector<double>& solution,
        const std::vector<double>& rhs,
        int nx,
        int ny,
        int nz,
        double hx,
        double hy,
        double hz
    ) const;

    // ------------------------------------------------------------
    // Suavizado Gauss-Seidel
    // ------------------------------------------------------------

    void smooth(
        std::vector<double>& solution,
        const std::vector<double>& rhs,
        int nx,
        int ny,
        int nz,
        double hx,
        double hy,
        double hz,
        int iterations
    ) const;

    // ------------------------------------------------------------
    // Residuo
    // ------------------------------------------------------------

    std::vector<double> calculateResidual(
        const std::vector<double>& solution,
        const std::vector<double>& rhs,
        int nx,
        int ny,
        int nz,
        double hx,
        double hy,
        double hz
    ) const;

    double calculateResidualNorm(
        const std::vector<double>& residual,
        int nx,
        int ny,
        int nz,
        double hx,
        double hy,
        double hz
    ) const;

    // ------------------------------------------------------------
    // Restricción
    // ------------------------------------------------------------

    std::vector<double> restrictResidual(
        const std::vector<double>& fineResidual,
        int fineNx,
        int fineNy,
        int fineNz,
        int coarseNx,
        int coarseNy,
        int coarseNz
    ) const;

    // ------------------------------------------------------------
    // Prolongación
    // ------------------------------------------------------------

    void prolongateAndAdd(
        const std::vector<double>& coarseCorrection,
        int coarseNx,
        int coarseNy,
        int coarseNz,
        std::vector<double>& fineSolution,
        int fineNx,
        int fineNy,
        int fineNz
    ) const;

    // ------------------------------------------------------------
    // Interpolación trilineal
    // ------------------------------------------------------------

    double sampleTrilinear(
        const std::vector<double>& field,
        int nx,
        int ny,
        int nz,
        double x,
        double y,
        double z
    ) const;
};