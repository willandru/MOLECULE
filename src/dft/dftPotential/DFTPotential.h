#pragma once

#include <vector>
#include <cstddef>

#include <glm/glm.hpp>

#include "DFTGrid.h"
#include "DFTDensity.h"

#include "ExternalPotential.h"
#include "HartreePotential.h"
#include "ExchangePotential.h"
#include "PotentialEnergy.h"
#include "NuclearForces.h"

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

    ExternalPotential externalPotential;

    HartreePotential hartreePotential;

    ExchangePotential exchangePotential;

    PotentialEnergy potentialEnergy;

    NuclearForces nuclearForces;
};