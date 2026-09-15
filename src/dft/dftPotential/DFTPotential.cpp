#include "DFTPotential.h"

#include <stdexcept>


// ================================================================
// CONSTRUCTORES
// ================================================================

DFTPotential::DFTPotential()
    :
    grid(nullptr),
    externalPotential(),
    hartreePotential(),
    exchangePotential(),
    potentialEnergy(),
    nuclearForces()
{
}


DFTPotential::DFTPotential(
    const DFTGrid& grid
)
    :
    grid(nullptr),
    externalPotential(),
    hartreePotential(),
    exchangePotential(),
    potentialEnergy(),
    nuclearForces()
{
    initialize(grid);
}


// ================================================================
// INICIALIZACIÓN
// ================================================================

void DFTPotential::initialize(
    const DFTGrid& grid
)
{
    if (grid.getPointCount() == 0)
    {
        throw std::invalid_argument(
            "DFTPotential: grid has no points."
        );
    }

    if (grid.getNx() < 3 ||
        grid.getNy() < 3 ||
        grid.getNz() < 3)
    {
        throw std::invalid_argument(
            "DFTPotential: grid dimensions must be at least 3 in every direction."
        );
    }

    if (!std::isfinite(grid.getSpacing()) ||
        grid.getSpacing() <= 0.0)
    {
        throw std::invalid_argument(
            "DFTPotential: grid spacing must be positive and finite."
        );
    }

    this->grid =
        &grid;

    externalPotential.initialize(
        grid
    );

    hartreePotential.initialize(
        grid
    );

    exchangePotential.initialize(
        grid
    );

    potentialEnergy.initialize(
        grid
    );

    nuclearForces.initialize(
        grid
    );
}


// ================================================================
// POTENCIAL EXTERNO
// ================================================================

void DFTPotential::calculateExternalPotential(
    const std::vector<int>& nuclearCharges,
    const std::vector<glm::dvec3>& nuclearPositions
)
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "DFTPotential: not initialized."
        );
    }

    externalPotential.calculate(
        nuclearCharges,
        nuclearPositions
    );
}


// ================================================================
// POTENCIAL DE HARTREE
// ================================================================

void DFTPotential::calculateHartreePotential(
    const DFTDensity& density,
    int iterations
)
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "DFTPotential: not initialized."
        );
    }

    hartreePotential.calculate(
        density,
        iterations
    );
}


// ================================================================
// POTENCIAL DE INTERCAMBIO
// ================================================================

void DFTPotential::calculateExchangePotential(
    const DFTDensity& density
)
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "DFTPotential: not initialized."
        );
    }

    exchangePotential.calculate(
        density
    );
}


// ================================================================
// ACCESO A POTENCIALES
// ================================================================

double DFTPotential::getExternal(
    std::size_t index
) const
{
    return externalPotential.get(
        index
    );
}


double DFTPotential::getHartree(
    std::size_t index
) const
{
    return hartreePotential.get(
        index
    );
}


double DFTPotential::getExchange(
    std::size_t index
) const
{
    return exchangePotential.get(
        index
    );
}


double DFTPotential::getTotalElectronic(
    std::size_t index
) const
{
    return
        externalPotential.get(index)
        +
        hartreePotential.get(index)
        +
        exchangePotential.get(index);
}


const std::vector<double>&
DFTPotential::getExternalPotential() const
{
    return externalPotential.getPotential();
}


const std::vector<double>&
DFTPotential::getHartreePotential() const
{
    return hartreePotential.getPotential();
}


const std::vector<double>&
DFTPotential::getExchangePotential() const
{
    return exchangePotential.getPotential();
}


// ================================================================
// DIAGNÓSTICO DEL SOLVER DE POISSON
// ================================================================

double DFTPotential::getPoissonResidual() const
{
    return hartreePotential.getResidual();
}


int DFTPotential::getPoissonCycles() const
{
    return hartreePotential.getCycles();
}


bool DFTPotential::hasPoissonConverged() const
{
    return hartreePotential.hasConverged();
}


// ================================================================
// ENERGÍA ELECTRÓN - NÚCLEO
// ================================================================

double DFTPotential::calculateElectronNuclearEnergy(
    const DFTDensity& density
) const
{
    return potentialEnergy.calculateElectronNuclearEnergy(
        density,
        externalPotential.getPotential()
    );
}


// ================================================================
// ENERGÍA DE HARTREE
// ================================================================

double DFTPotential::calculateHartreeEnergy(
    const DFTDensity& density
) const
{
    return potentialEnergy.calculateHartreeEnergy(
        density,
        hartreePotential.getPotential()
    );
}


// ================================================================
// ENERGÍA DE INTERCAMBIO
// ================================================================

double DFTPotential::calculateExchangeEnergy(
    const DFTDensity& density
) const
{
    return potentialEnergy.calculateExchangeEnergy(
        density
    );
}


// ================================================================
// REPULSIÓN NÚCLEO - NÚCLEO
// ================================================================

double DFTPotential::calculateNuclearRepulsionEnergy(
    const std::vector<int>& nuclearCharges,
    const std::vector<glm::dvec3>& nuclearPositions
) const
{
    return potentialEnergy.calculateNuclearRepulsionEnergy(
        nuclearCharges,
        nuclearPositions
    );
}


// ================================================================
// FUERZAS NUCLEARES
// ================================================================

std::vector<glm::dvec3>
DFTPotential::calculateNuclearForces(
    const DFTDensity& density,
    const std::vector<int>& nuclearCharges,
    const std::vector<glm::dvec3>& nuclearPositions
) const
{
    return nuclearForces.calculate(
        density,
        nuclearCharges,
        nuclearPositions,
        externalPotential.getSoftening()
    );
}


// ================================================================
// SOFTENING
// ================================================================

void DFTPotential::setSoftening(
    double value
)
{
    externalPotential.setSoftening(
        value
    );
}


double DFTPotential::getSoftening() const
{
    return externalPotential.getSoftening();
}


// ================================================================
// GRID
// ================================================================

const DFTGrid& DFTPotential::getGrid() const
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "DFTPotential: not initialized."
        );
    }

    return *grid;
}