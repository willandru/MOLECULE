#include "DFT.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>


namespace
{
    constexpr double MIN_ELECTRON_COUNT =
        1.0e-10;

    constexpr double MIN_DISTANCE =
        1.0e-12;
}


// ================================================================
// CONSTRUCTORES
// ================================================================

DFT::DFT()
    :
    grid(nullptr),
    electronNuclearEnergy(0.0),
    hartreeEnergy(0.0),
    exchangeEnergy(0.0),
    nuclearRepulsionEnergy(0.0),
    converged(false)
{
}


DFT::DFT(
    const DFTGrid& grid
)
    :
    grid(nullptr),
    electronNuclearEnergy(0.0),
    hartreeEnergy(0.0),
    exchangeEnergy(0.0),
    nuclearRepulsionEnergy(0.0),
    converged(false)
{
    initialize(grid);
}


// ================================================================
// INICIALIZACIÓN
// ================================================================

void DFT::initialize(
    const DFTGrid& grid
)
{
    if (grid.getPointCount() == 0)
    {
        throw std::invalid_argument(
            "DFT: grid has no points."
        );
    }

    this->grid =
        &grid;

    density.initialize(
        grid
    );

    potential.initialize(
        grid
    );

    kohnSham.initialize(
        grid
    );

    nuclearCharges.clear();

    nuclearPositions.clear();

    geometry.initialize(
        {},
        {}
    );

    electronNuclearEnergy =
        0.0;

    hartreeEnergy =
        0.0;

    exchangeEnergy =
        0.0;

    nuclearRepulsionEnergy =
        0.0;

    converged =
        false;
}


// ================================================================
// MOLÉCULA
// ================================================================

void DFT::setMolecule(
    const std::vector<int>& charges,
    const std::vector<glm::dvec3>& positions
)
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "DFT: not initialized."
        );
    }

    if (charges.empty())
    {
        throw std::invalid_argument(
            "DFT: molecule must contain at least one nucleus."
        );
    }

    if (charges.size() != positions.size())
    {
        throw std::invalid_argument(
            "DFT: charges and positions have different sizes."
        );
    }

    for (int charge : charges)
    {
        if (charge <= 0)
        {
            throw std::invalid_argument(
                "DFT: nuclear charge must be positive."
            );
        }
    }

    for (std::size_t i = 0;
         i < positions.size();
         ++i)
    {
        for (std::size_t j = i + 1;
             j < positions.size();
             ++j)
        {
            const glm::dvec3 difference =
                positions[i] -
                positions[j];

            const double distance =
                glm::length(
                    difference
                );

            if (distance <= MIN_DISTANCE)
            {
                throw std::invalid_argument(
                    "DFT: nuclei cannot occupy the same position."
                );
            }
        }
    }

    nuclearCharges =
        charges;

    nuclearPositions =
        positions;

    geometry.initialize(
        nuclearCharges,
        nuclearPositions
    );

    const int electronCount =
        calculateElectronCount();

    kohnSham.setElectronCount(
        electronCount
    );

    density.clear();

    initializeInitialDensity();

    converged =
        false;

    electronNuclearEnergy =
        0.0;

    hartreeEnergy =
        0.0;

    exchangeEnergy =
        0.0;

    nuclearRepulsionEnergy =
        0.0;
}


const std::vector<int>&
DFT::getNuclearCharges() const
{
    return nuclearCharges;
}


const std::vector<glm::dvec3>&
DFT::getNuclearPositions() const
{
    return nuclearPositions;
}


DFTGeometry&
DFT::getGeometry()
{
    return geometry;
}


const DFTGeometry&
DFT::getGeometry() const
{
    return geometry;
}


// ================================================================
// NÚMERO DE ELECTRONES
// ================================================================

int DFT::calculateElectronCount() const
{
    int electrons =
        0;

    for (int charge : nuclearCharges)
    {
        if (charge <= 0)
        {
            throw std::invalid_argument(
                "DFT: nuclear charge must be positive."
            );
        }

        electrons +=
            charge;
    }

    return electrons;
}


// ================================================================
// DENSIDAD INICIAL
// ================================================================

void DFT::initializeInitialDensity()
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "DFT: not initialized."
        );
    }

    if (nuclearCharges.empty())
    {
        throw std::runtime_error(
            "DFT: molecule has not been defined."
        );
    }

    const std::size_t count =
        grid->getPointCount();

    for (std::size_t i = 0;
         i < count;
         ++i)
    {
        const glm::dvec3 position =
            grid->getPosition(i);

        double value =
            0.0;

        for (std::size_t nucleus = 0;
             nucleus < nuclearCharges.size();
             ++nucleus)
        {
            const double Z =
                static_cast<double>(
                    nuclearCharges[nucleus]
                );

            const glm::dvec3 difference =
                position -
                nuclearPositions[nucleus];

            const double r2 =
                glm::dot(
                    difference,
                    difference
                );

            const double r =
                std::sqrt(
                    r2 +
                    0.01
                );

            value +=
                Z * Z * Z *
                std::exp(
                    -2.0 * Z * r
                );
        }

        density.set(
            i,
            value
        );
    }

    density.normalize(
        static_cast<double>(
            calculateElectronCount()
        )
    );
}


// ================================================================
// VALIDACIÓN DE DENSIDAD
// ================================================================

bool DFT::hasValidDensity() const
{
    if (grid == nullptr)
    {
        return false;
    }

    const double electrons =
        density.calculateElectronCount();

    const double expected =
        static_cast<double>(
            calculateElectronCount()
        );

    if (!std::isfinite(electrons))
    {
        return false;
    }

    return
        electrons > MIN_ELECTRON_COUNT &&
        std::abs(
            electrons - expected
        ) < 1.0e-6;
}


// ================================================================
// POTENCIALES
// ================================================================

void DFT::updatePotential()
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "DFT: not initialized."
        );
    }

    if (nuclearCharges.empty())
    {
        throw std::runtime_error(
            "DFT: molecule has not been defined."
        );
    }

    potential.calculateExternalPotential(
        nuclearCharges,
        nuclearPositions
    );

    potential.calculateHartreePotential(
        density
    );

    potential.calculateExchangePotential(
        density
    );
}


// ================================================================
// SCF
// ================================================================

bool DFT::solveSCF(
    int iterations,
    int orbitalIterations,
    double orbitalStep,
    double densityMixing,
    double densityTolerance
)
{
    std::cout
        << "[DFT] Entrando en solveSCF()"
        << std::endl;


    if (grid == nullptr)
    {
        throw std::runtime_error(
            "DFT: not initialized."
        );
    }


    if (nuclearCharges.empty())
    {
        throw std::runtime_error(
            "DFT: molecule has not been defined."
        );
    }


    if (iterations <= 0)
    {
        throw std::invalid_argument(
            "DFT: SCF iterations must be positive."
        );
    }


    if (orbitalIterations <= 0)
    {
        throw std::invalid_argument(
            "DFT: orbital iterations must be positive."
        );
    }


    if (!std::isfinite(orbitalStep) ||
        orbitalStep <= 0.0)
    {
        throw std::invalid_argument(
            "DFT: orbital step must be positive and finite."
        );
    }


    if (!std::isfinite(densityMixing) ||
        densityMixing <= 0.0 ||
        densityMixing > 1.0)
    {
        throw std::invalid_argument(
            "DFT: density mixing must be in (0,1]."
        );
    }


    if (!std::isfinite(densityTolerance) ||
        densityTolerance <= 0.0)
    {
        throw std::invalid_argument(
            "DFT: density tolerance must be positive and finite."
        );
    }


    // ============================================================
    // DENSIDAD INICIAL
    // ============================================================

    std::cout
        << "[DFT] Comprobando densidad..."
        << std::endl;


    if (!hasValidDensity())
    {
        std::cout
            << "[DFT] Densidad no valida."
            << std::endl;


        density.clear();


        std::cout
            << "[DFT] Construyendo densidad inicial..."
            << std::endl;


        initializeInitialDensity();


        std::cout
            << "[DFT] Densidad inicial construida."
            << std::endl;
    }
    else
    {
        std::cout
            << "[DFT] Densidad existente reutilizada."
            << std::endl;
    }


    converged =
        false;


    // ============================================================
    // CICLO SCF
    // ============================================================

    for (int iteration = 0;
         iteration < iterations;
         ++iteration)
    {
        std::cout
            << "\n[DFT] SCF iteration "
            << iteration + 1
            << " / "
            << iterations
            << std::endl;


        // --------------------------------------------------------
        // DENSIDAD ANTERIOR
        // --------------------------------------------------------

        std::cout
            << "[DFT] Copiando densidad anterior..."
            << std::endl;


        const DFTDensity oldDensity =
            density;


        std::cout
            << "[DFT] Densidad anterior copiada."
            << std::endl;


        // --------------------------------------------------------
        // POTENCIALES
        // --------------------------------------------------------

        std::cout
            << "[DFT] Calculando potenciales..."
            << std::endl;


        updatePotential();


        std::cout
            << "[DFT] Potenciales calculados."
            << std::endl;


        // --------------------------------------------------------
        // KOHN-SHAM
        // --------------------------------------------------------

        std::cout
            << "[DFT] Entrando en Kohn-Sham..."
            << std::endl;


        kohnSham.solve(
            potential,
            orbitalIterations,
            orbitalStep
        );


        std::cout
            << "[DFT] Kohn-Sham terminado."
            << std::endl;


        // --------------------------------------------------------
        // NUEVA DENSIDAD
        // --------------------------------------------------------

        std::cout
            << "[DFT] Construyendo densidad desde orbitales..."
            << std::endl;


        DFTDensity calculatedDensity(
            *grid
        );


        std::cout
            << "[DFT] Objeto calculatedDensity creado."
            << std::endl;


        kohnSham.calculateDensity(
            calculatedDensity
        );


        std::cout
            << "[DFT] Densidad KS construida."
            << std::endl;


        // --------------------------------------------------------
        // NORMALIZACIÓN
        // --------------------------------------------------------

        std::cout
            << "[DFT] Normalizando densidad KS..."
            << std::endl;


        calculatedDensity.normalize(
            static_cast<double>(
                calculateElectronCount()
            )
        );


        std::cout
            << "[DFT] Densidad KS normalizada."
            << std::endl;


        // --------------------------------------------------------
        // MEZCLA
        // --------------------------------------------------------

        std::cout
            << "[DFT] Mezclando densidades..."
            << std::endl;


        density.mix(
            calculatedDensity,
            densityMixing
        );


        std::cout
            << "[DFT] Mezcla terminada."
            << std::endl;


        // --------------------------------------------------------
        // NORMALIZACIÓN FINAL
        // --------------------------------------------------------

        std::cout
            << "[DFT] Normalizando densidad final..."
            << std::endl;


        density.normalize(
            static_cast<double>(
                calculateElectronCount()
            )
        );


        std::cout
            << "[DFT] Normalización final terminada."
            << std::endl;


        // --------------------------------------------------------
        // CONVERGENCIA
        // --------------------------------------------------------

        std::cout
            << "[DFT] Calculando diferencia de densidad..."
            << std::endl;


        const double densityDifference =
            density.calculateDifferenceNorm(
                oldDensity
            );


        std::cout
            << "[DFT] Density difference = "
            << densityDifference
            << std::endl;


        if (!std::isfinite(
                densityDifference
            ))
        {
            throw std::runtime_error(
                "DFT: non-finite SCF density difference."
            );
        }


        if (densityDifference <
                densityTolerance &&
            kohnSham.hasConverged())
        {
            converged =
                true;


            std::cout
                << "[DFT] SCF CONVERGED."
                << std::endl;


            break;
        }


        std::cout
            << "[DFT] SCF aun no convergido."
            << std::endl;
    }


    // ============================================================
    // POTENCIAL FINAL
    // ============================================================

    std::cout
        << "\n[DFT] Calculando potencial final..."
        << std::endl;


    updatePotential();


    std::cout
        << "[DFT] Potencial final calculado."
        << std::endl;


    // ============================================================
    // ENERGÍA FINAL
    // ============================================================

    std::cout
        << "[DFT] Calculando energias..."
        << std::endl;


    updateEnergy();


    std::cout
        << "[DFT] Energias calculadas."
        << std::endl;


    std::cout
        << "[DFT] solveSCF() terminado."
        << std::endl;


    return converged;
}


// ================================================================
// ENERGÍA CINÉTICA
// ================================================================

double DFT::calculateKineticEnergy() const
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "DFT: not initialized."
        );
    }

    const int nx =
        grid->getNx();

    const int ny =
        grid->getNy();

    const int nz =
        grid->getNz();

    const double h =
        grid->getSpacing();

    const double h2 =
        h * h;

    const double dV =
        grid->getVolumeElement();

    double kineticEnergy =
        0.0;

    const int orbitalCount =
        kohnSham.getOrbitalCount();


    for (int orbital = 0;
         orbital < orbitalCount;
         ++orbital)
    {
        const double occupation =
            kohnSham.getOccupation(
                orbital
            );


        if (occupation <= 0.0)
        {
            continue;
        }


        const std::vector<double>& psi =
            kohnSham.getOrbital(
                orbital
            );


        double orbitalKinetic =
            0.0;


        for (int z = 1;
             z < nz - 1;
             ++z)
        {
            for (int y = 1;
                 y < ny - 1;
                 ++y)
            {
                for (int x = 1;
                     x < nx - 1;
                     ++x)
                {
                    const std::size_t index =
                        grid->getIndex(
                            x,
                            y,
                            z
                        );


                    const double neighbors =
                        psi[
                            grid->getIndex(
                                x + 1,
                                y,
                                z
                            )
                        ]
                        +
                        psi[
                            grid->getIndex(
                                x - 1,
                                y,
                                z
                            )
                        ]
                        +
                        psi[
                            grid->getIndex(
                                x,
                                y + 1,
                                z
                            )
                        ]
                        +
                        psi[
                            grid->getIndex(
                                x,
                                y - 1,
                                z
                            )
                        ]
                        +
                        psi[
                            grid->getIndex(
                                x,
                                y,
                                z + 1
                            )
                        ]
                        +
                        psi[
                            grid->getIndex(
                                x,
                                y,
                                z - 1
                            )
                        ];


                    const double laplacian =
                        (
                            neighbors -
                            6.0 *
                            psi[index]
                        )
                        /
                        h2;


                    orbitalKinetic +=
                        psi[index] *
                        (
                            -0.5 *
                            laplacian
                        );
                }
            }
        }


        orbitalKinetic *=
            dV;


        kineticEnergy +=
            occupation *
            orbitalKinetic;
    }


    return kineticEnergy;
}


// ================================================================
// ENERGÍA
// ================================================================

void DFT::updateEnergy()
{
    electronNuclearEnergy =
        potential.calculateElectronNuclearEnergy(
            density
        );

    hartreeEnergy =
        potential.calculateHartreeEnergy(
            density
        );

    exchangeEnergy =
        potential.calculateExchangeEnergy(
            density
        );

    nuclearRepulsionEnergy =
        potential.calculateNuclearRepulsionEnergy(
            nuclearCharges,
            nuclearPositions
        );
}


double DFT::calculateTotalEnergy() const
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "DFT: not initialized."
        );
    }

    const double kineticEnergy =
        calculateKineticEnergy();

    return
        kineticEnergy
        +
        electronNuclearEnergy
        +
        hartreeEnergy
        +
        exchangeEnergy
        +
        nuclearRepulsionEnergy;
}


double DFT::getElectronNuclearEnergy() const
{
    return electronNuclearEnergy;
}


double DFT::getHartreeEnergy() const
{
    return hartreeEnergy;
}


double DFT::getExchangeEnergy() const
{
    return exchangeEnergy;
}


double DFT::getNuclearRepulsionEnergy() const
{
    return nuclearRepulsionEnergy;
}


// ================================================================
// FUERZAS
// ================================================================

std::vector<glm::dvec3>
DFT::calculateForces() const
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "DFT: not initialized."
        );
    }

    if (nuclearCharges.empty())
    {
        throw std::runtime_error(
            "DFT: molecule has not been defined."
        );
    }

    return
        potential.calculateNuclearForces(
            density,
            nuclearCharges,
            nuclearPositions
        );
}


// ================================================================
// OPTIMIZACIÓN GEOMÉTRICA
// ================================================================

// ================================================================
// OPTIMIZACIÓN GEOMÉTRICA
// ================================================================

bool DFT::optimizeGeometry(
    int geometryIterations,
    int scfIterations,
    int orbitalIterations,
    double orbitalStep,
    double densityMixing,
    double densityTolerance,
    double geometryStep,
    double forceTolerance
)
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "DFT: not initialized."
        );
    }

    if (nuclearCharges.empty())
    {
        throw std::runtime_error(
            "DFT: molecule must be defined."
        );
    }

    if (geometryIterations <= 0)
    {
        throw std::invalid_argument(
            "DFT: geometry iterations must be positive."
        );
    }

    if (!std::isfinite(geometryStep) ||
        geometryStep <= 0.0)
    {
        throw std::invalid_argument(
            "DFT: geometry step must be positive and finite."
        );
    }

    if (!std::isfinite(forceTolerance) ||
        forceTolerance <= 0.0)
    {
        throw std::invalid_argument(
            "DFT: force tolerance must be positive and finite."
        );
    }


    for (int iteration = 0;
         iteration < geometryIterations;
         ++iteration)
    {
        std::cout
            << "\n[DFT] Geometry iteration "
            << iteration + 1
            << " / "
            << geometryIterations
            << std::endl;


        // ========================================================
        // SCF ELECTRÓNICO
        // ========================================================

        const bool electronicConverged =
            solveSCF(
                scfIterations,
                orbitalIterations,
                orbitalStep,
                densityMixing,
                densityTolerance
            );


        // ========================================================
        // SI EL SCF NO CONVERGIÓ, DETENER TODO
        // ========================================================

        if (!electronicConverged)
        {
            std::cout
                << "[DFT] ERROR: el SCF no convergio."
                << std::endl;

            std::cout
                << "[DFT] La optimizacion geometrica se detiene."
                << std::endl;

            return false;
        }


        // ========================================================
        // CALCULAR FUERZAS
        // ========================================================

        std::cout
            << "[DFT] Calculando fuerzas..."
            << std::endl;


        const std::vector<glm::dvec3> forces =
            calculateForces();


        geometry.setForces(
            forces
        );


        // ========================================================
        // COMPROBAR CONVERGENCIA GEOMÉTRICA
        // ========================================================

        const double maximumForce =
            geometry.getMaximumForce();


        std::cout
            << "[DFT] Maximum force = "
            << maximumForce
            << " Ha/Bohr"
            << std::endl;


        if (geometry.isOptimized(
                forceTolerance
            ))
        {
            std::cout
                << "[DFT] GEOMETRIA CONVERGIDA."
                << std::endl;

            return true;
        }


        // ========================================================
        // MOVER NÚCLEOS
        // ========================================================

        std::cout
            << "[DFT] Moviendo nucleos..."
            << std::endl;


        geometry.moveNuclei(
            geometryStep
        );


        // ========================================================
        // ACTUALIZAR POSICIONES DE DFT
        // ========================================================

        nuclearPositions =
            geometry.getNuclearPositions();


        converged =
            false;
    }


    // ============================================================
    // SE AGOTARON LAS ITERACIONES GEOMÉTRICAS
    // ============================================================

    const double maximumForce =
        geometry.getMaximumForce();


    std::cout
        << "[DFT] Se alcanzo el maximo de iteraciones geometricas."
        << std::endl;


    std::cout
        << "[DFT] Maximum force = "
        << maximumForce
        << " Ha/Bohr"
        << std::endl;


    return
        geometry.isOptimized(
            forceTolerance
        );
}

// ================================================================
// RESULTADOS
// ================================================================

const DFTDensity&
DFT::getDensity() const
{
    return density;
}


const DFTKohnSham&
DFT::getKohnSham() const
{
    return kohnSham;
}


const DFTPotential&
DFT::getPotential() const
{
    return potential;
}


bool DFT::hasConverged() const
{
    return converged;
}


// ================================================================
// CONFIGURACIÓN
// ================================================================

void DFT::setPotentialSoftening(
    double value
)
{
    potential.setSoftening(
        value
    );
}


// ================================================================
// GRID
// ================================================================

const DFTGrid&
DFT::getGrid() const
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "DFT: not initialized."
        );
    }

    return *grid;
}