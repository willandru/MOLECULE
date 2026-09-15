#pragma once

#include <vector>

#include <glm/glm.hpp>

#include "DFTGrid.h"
#include "DFTDensity.h"
#include "DFTPotential.h"
#include "DFTKohnSham.h"
#include "DFTGeometry.h"


class DFT
{
public:

    DFT();

    explicit DFT(
        const DFTGrid& grid
    );

    void initialize(
        const DFTGrid& grid
    );


    // ============================================================
    // MOLÉCULA
    // ============================================================

    void setMolecule(
        const std::vector<int>& nuclearCharges,
        const std::vector<glm::dvec3>& nuclearPositions
    );

    const std::vector<int>& getNuclearCharges() const;

    const std::vector<glm::dvec3>& getNuclearPositions() const;

    DFTGeometry& getGeometry();

    const DFTGeometry& getGeometry() const;


    // ============================================================
    // SCF
    // ============================================================

    bool solveSCF(
        int iterations = 100,
        int orbitalIterations = 200,
        double orbitalStep = 0.001,
        double densityMixing = 0.25,
        double densityTolerance = 1.0e-6
    );


    // ============================================================
    // ENERGÍA
    // ============================================================

    double calculateTotalEnergy() const;

    double getElectronNuclearEnergy() const;

    double getHartreeEnergy() const;

    double getExchangeEnergy() const;

    double getNuclearRepulsionEnergy() const;


    // ============================================================
    // FUERZAS
    // ============================================================

    std::vector<glm::dvec3> calculateForces() const;


    // ============================================================
    // OPTIMIZACIÓN GEOMÉTRICA
    // ============================================================

    bool optimizeGeometry(
        int geometryIterations = 50,
        int scfIterations = 100,
        int orbitalIterations = 200,
        double orbitalStep = 0.001,
        double densityMixing = 0.25,
        double densityTolerance = 1.0e-6,
        double geometryStep = 0.05,
        double forceTolerance = 1.0e-4
    );


    // ============================================================
    // RESULTADOS ELECTRÓNICOS
    // ============================================================

    const DFTDensity& getDensity() const;

    const DFTKohnSham& getKohnSham() const;

    const DFTPotential& getPotential() const;

    bool hasConverged() const;


    // ============================================================
    // CONFIGURACIÓN
    // ============================================================

    void setPotentialSoftening(
        double value
    );

    const DFTGrid& getGrid() const;


private:

    const DFTGrid* grid;

    DFTDensity density;

    DFTPotential potential;

    DFTKohnSham kohnSham;

    DFTGeometry geometry;

    std::vector<int> nuclearCharges;

    std::vector<glm::dvec3> nuclearPositions;

    double electronNuclearEnergy;

    double hartreeEnergy;

    double exchangeEnergy;

    double nuclearRepulsionEnergy;

    bool converged;


    // ============================================================
    // MÉTODOS INTERNOS
    // ============================================================

    void initializeInitialDensity();

    void updatePotential();

    void updateEnergy();

    double calculateKineticEnergy() const;

    int calculateElectronCount() const;

    bool hasValidDensity() const;
};