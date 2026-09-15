#pragma once

#include <vector>
#include <cstddef>

#include <glm/glm.hpp>

#include "DFTGrid.h"
#include "DFTDensity.h"

class DFTPotential
{
public:

    DFTPotential();

    explicit DFTPotential(
        const DFTGrid& grid
    );

    void initialize(
        const DFTGrid& grid
    );

    // ------------------------------------------------------------
    // Potencial externo electrón - núcleo
    // ------------------------------------------------------------

    void calculateExternalPotential(
        const std::vector<int>& nuclearCharges,
        const std::vector<glm::dvec3>& nuclearPositions
    );

    // ------------------------------------------------------------
    // Potencial de Hartree electrón - electrón
    // ------------------------------------------------------------

    void calculateHartreePotential(
        const DFTDensity& density,
        int iterations = 20
    );

    // ------------------------------------------------------------
    // Potencial de intercambio LDA
    // ------------------------------------------------------------

    void calculateExchangePotential(
        const DFTDensity& density
    );

    // ------------------------------------------------------------
    // Acceso
    // ------------------------------------------------------------

    double getExternal(
        std::size_t index
    ) const;

    double getHartree(
        std::size_t index
    ) const;

    double getExchange(
        std::size_t index
    ) const;

    double getTotalElectronic(
        std::size_t index
    ) const;

    const std::vector<double>& getExternalPotential() const;
    const std::vector<double>& getHartreePotential() const;
    const std::vector<double>& getExchangePotential() const;

    // ------------------------------------------------------------
    // Diagnóstico del solver de Poisson
    // ------------------------------------------------------------

    double getPoissonResidual() const;

    int getPoissonCycles() const;

    bool hasPoissonConverged() const;

    // ------------------------------------------------------------
    // Energías
    // ------------------------------------------------------------

    double calculateElectronNuclearEnergy(
        const DFTDensity& density
    ) const;

    double calculateHartreeEnergy(
        const DFTDensity& density
    ) const;

    double calculateExchangeEnergy(
        const DFTDensity& density
    ) const;

    double calculateNuclearRepulsionEnergy(
        const std::vector<int>& nuclearCharges,
        const std::vector<glm::dvec3>& nuclearPositions
    ) const;

    // ------------------------------------------------------------
    // Fuerzas nucleares
    // ------------------------------------------------------------

    std::vector<glm::dvec3> calculateNuclearForces(
        const DFTDensity& density,
        const std::vector<int>& nuclearCharges,
        const std::vector<glm::dvec3>& nuclearPositions
    ) const;

    // ------------------------------------------------------------
    // Configuración
    // ------------------------------------------------------------

    void setSoftening(
        double value
    );

    double getSoftening() const;

    const DFTGrid& getGrid() const;

private:

    const DFTGrid* grid;

    std::vector<double> externalPotential;
    std::vector<double> hartreePotential;
    std::vector<double> exchangePotential;

    double softening;

    // ------------------------------------------------------------
    // Diagnóstico de Poisson
    // ------------------------------------------------------------

    double poissonResidual;

    int poissonCycles;

    bool poissonConverged;

    // ------------------------------------------------------------
    // Poisson / Multigrid
    // ------------------------------------------------------------

    void solvePoisson(
        const DFTDensity& density,
        std::vector<double>& potential,
        int maxCycles
    );

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

    std::vector<double> restrictResidual(
        const std::vector<double>& fineResidual,
        int fineNx,
        int fineNy,
        int fineNz,
        int coarseNx,
        int coarseNy,
        int coarseNz
    ) const;

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

    double sampleTrilinear(
        const std::vector<double>& field,
        int nx,
        int ny,
        int nz,
        double x,
        double y,
        double z
    ) const;

    double calculateLaplacian(
        const std::vector<double>& potential,
        int x,
        int y,
        int z
    ) const;

    double calculateDensityCharge(
        const DFTDensity& density
    ) const;

    // ------------------------------------------------------------
    // Condiciones de frontera
    // ------------------------------------------------------------

    void applyBoundaryConditions(
        std::vector<double>& potential,
        const DFTDensity& density
    ) const;
};