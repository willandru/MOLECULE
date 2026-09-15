#include "DFTHydrogen.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>


namespace
{
    constexpr double PI =
        3.1415926535897932384626433832795;

    constexpr double FOUR_PI =
        4.0 * PI;


    /*
        Unidades atómicas:

            hbar = 1
            me   = 1
            e    = 1
            a0   = 1

        Para hidrógeno:

            V_nuclear(r) = -1/r

        Intercambio LDA:

            V_x(r) = -(3/pi)^(1/3) rho(r)^(1/3)
    */
}


// ============================================================
// CONSTRUCTORES
// ============================================================

DFTHydrogen::DFTHydrogen()
    : parameters(),
      result()
{
}


DFTHydrogen::DFTHydrogen(
    const Parameters& parameters
)
    : parameters(parameters),
      result()
{
    if (this->parameters.gridPoints < 16)
    {
        throw std::invalid_argument(
            "DFTHydrogen: gridPoints debe ser >= 16."
        );
    }


    if (this->parameters.rMax <= 0.0)
    {
        throw std::invalid_argument(
            "DFTHydrogen: rMax debe ser > 0."
        );
    }


    if (this->parameters.mixing <= 0.0 ||
        this->parameters.mixing > 1.0)
    {
        throw std::invalid_argument(
            "DFTHydrogen: mixing debe estar en (0,1]."
        );
    }


    if (this->parameters.convergence <= 0.0)
    {
        throw std::invalid_argument(
            "DFTHydrogen: convergence debe ser > 0."
        );
    }


    if (this->parameters.maxSCFIterations == 0)
    {
        throw std::invalid_argument(
            "DFTHydrogen: maxSCFIterations debe ser > 0."
        );
    }
}


// ============================================================
// SOLVE
// ============================================================

DFTHydrogen::Result DFTHydrogen::solve()
{
    result = Result();


    std::cout
        << "\n"
        << "============================================================\n"
        << "                 DFT - ATOMO DE HIDROGENO\n"
        << "============================================================\n"
        << "\n";


    // ========================================================
    // PARAMETERS
    // ========================================================

    std::cout
        << "[DFT] Parametros\n"
        << "      Grid points       : "
        << parameters.gridPoints
        << "\n"
        << "      Radio maximo      : "
        << parameters.rMax
        << " Bohr\n"
        << "      Mixing            : "
        << parameters.mixing
        << "\n"
        << "      Convergencia      : "
        << std::scientific
        << parameters.convergence
        << "\n"
        << "      Max SCF iterations: "
        << parameters.maxSCFIterations
        << "\n"
        << std::defaultfloat
        << "\n";


    // ========================================================
    // RADIAL GRID
    // ========================================================

    std::cout
        << "[DFT] Construyendo malla radial..."
        << std::endl;


    buildRadialGrid();


    const std::size_t N =
        result.r.size();


    if (N < 16)
    {
        throw std::runtime_error(
            "DFTHydrogen: la malla radial es demasiado pequena."
        );
    }


    std::cout
        << "[DFT] Malla radial creada\n"
        << "      Puntos            : "
        << N
        << "\n"
        << "      r inicial         : "
        << result.r.front()
        << " Bohr\n"
        << "      r final           : "
        << result.r.back()
        << " Bohr\n"
        << "      dr                : "
        << result.r[1] - result.r[0]
        << " Bohr\n"
        << "\n";


    // ========================================================
    // INITIAL ORBITAL
    // ========================================================

    std::cout
        << "[DFT] Inicializando orbital 1s..."
        << std::endl;


    initializeOrbital();


    std::cout
        << "[DFT] Orbital inicial normalizado."
        << "\n"
        << "\n";


    // ========================================================
    // INITIAL DENSITY
    // ========================================================

    std::vector<double> density(
        N,
        0.0
    );


    calculateDensity(
        result.orbital,
        density
    );


    result.density =
        density;


    double maximumInitialDensity =
        0.0;


    for (double value : density)
    {
        maximumInitialDensity =
            std::max(
                maximumInitialDensity,
                value
            );
    }


    std::cout
        << "[DFT] Densidad inicial calculada\n"
        << "      Densidad maxima  : "
        << std::scientific
        << maximumInitialDensity
        << " e-/Bohr^3\n"
        << std::defaultfloat
        << "\n";


    // ========================================================
    // SCF
    // ========================================================

    std::vector<double> oldDensity =
        density;


    std::vector<double> newDensity(
        N,
        0.0
    );


    std::vector<double> mixedDensity(
        N,
        0.0
    );


    std::vector<double> hartreePotential(
        N,
        0.0
    );


    std::vector<double> exchangePotential(
        N,
        0.0
    );


    std::vector<double> newOrbital =
        result.orbital;


    std::cout
        << "------------------------------------------------------------\n"
        << " SCF\n"
        << "------------------------------------------------------------\n";


    std::cout
        << std::setw(8)
        << "Iter"
        << std::setw(20)
        << "Eigenvalue"
        << std::setw(20)
        << "Density error"
        << "\n";


    std::cout
        << "------------------------------------------------------------\n";


    for (std::size_t iteration = 0;
         iteration < parameters.maxSCFIterations;
         ++iteration)
    {
        // ====================================================
        // HARTREE
        // ====================================================

        calculateHartreePotential(
            oldDensity,
            hartreePotential
        );


        // ====================================================
        // EXCHANGE LDA
        // ====================================================

        calculateExchangePotential(
            oldDensity,
            exchangePotential
        );


        // ====================================================
        // KOHN-SHAM POTENTIAL
        // ====================================================

        calculateEffectivePotential(
            hartreePotential,
            exchangePotential,
            result.effectivePotential
        );


        // ====================================================
        // SOLVE KOHN-SHAM
        // ====================================================

        const double eigenvalue =
            solveKohnShamOrbital(
                result.effectivePotential,
                newOrbital
            );


        // ====================================================
        // NEW DENSITY
        // ====================================================

        calculateDensity(
            newOrbital,
            newDensity
        );


        // ====================================================
        // DENSITY MIXING
        // ====================================================

        mixDensity(
            oldDensity,
            newDensity,
            mixedDensity
        );


        // ====================================================
        // CONVERGENCE ERROR
        // ====================================================

        const double densityDifference =
            calculateDensityDifference(
                oldDensity,
                mixedDensity
            );


        std::cout
            << std::setw(8)
            << iteration + 1
            << std::setw(20)
            << std::scientific
            << std::setprecision(10)
            << eigenvalue
            << std::setw(20)
            << densityDifference
            << std::defaultfloat
            << "\n";


        // ====================================================
        // SAVE CURRENT STATE
        // ====================================================

        oldDensity =
            mixedDensity;


        result.orbital =
            newOrbital;


        result.density =
            oldDensity;


        result.eigenvalue =
            eigenvalue;


        result.iterations =
            iteration + 1;


        // ====================================================
        // CONVERGENCE
        // ====================================================

        if (densityDifference <
            parameters.convergence)
        {
            result.converged =
                true;


            std::cout
                << "------------------------------------------------------------\n"
                << "[DFT] SCF CONVERGIO\n"
                << "      Iteraciones: "
                << result.iterations
                << "\n"
                << "      Error      : "
                << std::scientific
                << densityDifference
                << "\n"
                << std::defaultfloat
                << "------------------------------------------------------------\n"
                << "\n";


            break;
        }
    }


    if (!result.converged)
    {
        std::cout
            << "------------------------------------------------------------\n"
            << "[DFT] ADVERTENCIA: SCF NO CONVERGIO\n"
            << "      Iteraciones realizadas: "
            << result.iterations
            << "\n"
            << "------------------------------------------------------------\n"
            << "\n";
    }


    // ========================================================
    // FINAL POTENTIALS
    // ========================================================

    calculateHartreePotential(
        result.density,
        hartreePotential
    );


    calculateExchangePotential(
        result.density,
        exchangePotential
    );


    calculateEffectivePotential(
        hartreePotential,
        exchangePotential,
        result.effectivePotential
    );


    // ========================================================
    // ENERGIES
    // ========================================================

    result.kineticEnergy =
        calculateKineticEnergy(
            result.orbital
        );


    result.nuclearAttractionEnergy =
        calculateNuclearAttractionEnergy(
            result.density
        );


    result.hartreeEnergy =
        calculateHartreeEnergy(
            result.density,
            hartreePotential
        );


    result.exchangeEnergy =
        calculateExchangeEnergy(
            result.density
        );


    result.totalEnergy =
        calculateTotalEnergy(
            result.kineticEnergy,
            result.nuclearAttractionEnergy,
            result.hartreeEnergy,
            result.exchangeEnergy
        );


    // ========================================================
    // FINAL DENSITY INFORMATION
    // ========================================================

    double maximumDensity =
        0.0;


    std::size_t maximumDensityIndex =
        0;


    for (std::size_t i = 0;
         i < result.density.size();
         ++i)
    {
        if (result.density[i] >
            maximumDensity)
        {
            maximumDensity =
                result.density[i];

            maximumDensityIndex =
                i;
        }
    }


    const double densityMaximumRadius =
        result.r[maximumDensityIndex];


    // ========================================================
    // ENERGY REPORT
    // ========================================================

    std::cout
        << "============================================================\n"
        << "                 RESULTADO DFT\n"
        << "============================================================\n"
        << "\n";


    std::cout
        << "Convergencia\n"
        << "------------\n"
        << "  Convergio           : "
        << (result.converged ? "SI" : "NO")
        << "\n"
        << "  Iteraciones SCF     : "
        << result.iterations
        << "\n"
        << "\n";


    std::cout
        << "Orbital 1s\n"
        << "----------\n"
        << "  Eigenvalue KS       : "
        << std::scientific
        << std::setprecision(12)
        << result.eigenvalue
        << " Ha\n"
        << "\n";


    std::cout
        << "Densidad electronica\n"
        << "--------------------\n"
        << "  Max rho             : "
        << maximumDensity
        << " e-/Bohr^3\n"
        << "  r(max rho)          : "
        << densityMaximumRadius
        << " Bohr\n"
        << "\n";


    std::cout
        << "Energias\n"
        << "--------\n"
        << "  Energia cinetica    : "
        << result.kineticEnergy
        << " Ha\n"
        << "  Nuclear-electron    : "
        << result.nuclearAttractionEnergy
        << " Ha\n"
        << "  Hartree             : "
        << result.hartreeEnergy
        << " Ha\n"
        << "  Exchange LDA        : "
        << result.exchangeEnergy
        << " Ha\n"
        << "  -------------------------------\n"
        << "  Energia total       : "
        << result.totalEnergy
        << " Ha\n"
        << "\n"
        << std::defaultfloat;


    // ========================================================
    // DFT -> GEOMETRY INFORMATION
    // ========================================================

    const double isovalue =
        0.02;


    const double targetDensity =
        maximumDensity *
        isovalue;


    double isosurfaceRadius =
        result.r.back();


    for (std::size_t i = 1;
         i < result.r.size();
         ++i)
    {
        const double d0 =
            result.density[i - 1];


        const double d1 =
            result.density[i];


        if (d0 >= targetDensity &&
            d1 <= targetDensity)
        {
            const double r0 =
                result.r[i - 1];


            const double r1 =
                result.r[i];


            const double denominator =
                d1 - d0;


            if (std::abs(denominator) <
                1.0e-15)
            {
                isosurfaceRadius =
                    r1;
            }
            else
            {
                const double t =
                    (targetDensity - d0) /
                    denominator;


                isosurfaceRadius =
                    r0 +
                    t *
                    (r1 - r0);
            }


            break;
        }
    }


    constexpr double BOHR_TO_ANGSTROM =
        0.529177210903;


    const double isosurfaceRadiusAngstrom =
        isosurfaceRadius *
        BOHR_TO_ANGSTROM;


    std::cout
        << "Isosuperficie para visualizacion\n"
        << "---------------------------------\n"
        << "  Isovalue            : "
        << isovalue
        << " * rho_max\n"
        << "  Target density      : "
        << std::scientific
        << targetDensity
        << " e-/Bohr^3\n"
        << "  Radio isosuperficie : "
        << isosurfaceRadius
        << " Bohr\n"
        << "  Radio geometrico    : "
        << isosurfaceRadiusAngstrom
        << " Angstrom\n"
        << "\n";


    std::cout
        << "============================================================\n"
        << "                 FIN DEL CALCULO DFT\n"
        << "============================================================\n"
        << std::endl;


    return result;
}


// ============================================================
// RESULTADO
// ============================================================

const DFTHydrogen::Result&
DFTHydrogen::getResult() const
{
    return result;
}


// ============================================================
// MALLA RADIAL
// ============================================================

void DFTHydrogen::buildRadialGrid()
{
    const std::size_t N =
        parameters.gridPoints;


    result.r.resize(N);


    /*
        No usamos r = 0 exactamente porque existen términos:

            1/r

        y:

            rho = u² / (4 pi r²)
    */

    const double dr =
        parameters.rMax /
        static_cast<double>(
            N - 1
        );


    for (std::size_t i = 0;
         i < N;
         ++i)
    {
        result.r[i] =
            static_cast<double>(
                i + 1
            ) *
            dr;
    }
}


// ============================================================
// ORBITAL INICIAL
// ============================================================

void DFTHydrogen::initializeOrbital()
{
    const std::size_t N =
        result.r.size();


    result.orbital.resize(N);


    /*
        Estado 1s:

            psi_1s(r) ∝ exp(-r)

        Orbital radial reducido:

            u(r) = r exp(-r)
    */

    for (std::size_t i = 0;
         i < N;
         ++i)
    {
        const double r =
            result.r[i];


        result.orbital[i] =
            r *
            std::exp(-r);
    }


    normalizeOrbital(
        result.orbital
    );
}


// ============================================================
// NORMALIZACIÓN
// ============================================================

void DFTHydrogen::normalizeOrbital(
    std::vector<double>& orbital
) const
{
    const double n =
        norm(orbital);


    if (n <=
        std::numeric_limits<double>::epsilon())
    {
        throw std::runtime_error(
            "DFTHydrogen: no se puede normalizar "
            "un orbital nulo."
        );
    }


    for (double& value : orbital)
    {
        value /= n;
    }
}


// ============================================================
// DENSIDAD
// ============================================================

void DFTHydrogen::calculateDensity(
    const std::vector<double>& orbital,
    std::vector<double>& density
) const
{
    const std::size_t N =
        result.r.size();


    density.resize(N);


    /*
        psi(r) = u(r) / (sqrt(4 pi) r)

        rho(r) = |psi(r)|²
               = u(r)² / (4 pi r²)
    */

    for (std::size_t i = 0;
         i < N;
         ++i)
    {
        const double r =
            result.r[i];


        if (r <= 0.0)
        {
            density[i] =
                0.0;

            continue;
        }


        density[i] =
            orbital[i] *
            orbital[i] /
            (FOUR_PI * r * r);
    }
}


// ============================================================
// POTENCIAL DE HARTREE
// ============================================================

void DFTHydrogen::calculateHartreePotential(
    const std::vector<double>& density,
    std::vector<double>& potential
) const
{
    const std::size_t N =
        result.r.size();


    potential.assign(
        N,
        0.0
    );


    std::vector<double> innerCharge(
        N,
        0.0
    );


    std::vector<double> outerContribution(
        N,
        0.0
    );


    for (std::size_t i = 1;
         i < N;
         ++i)
    {
        const double r0 =
            result.r[i - 1];


        const double r1 =
            result.r[i];


        const double f0 =
            FOUR_PI *
            density[i - 1] *
            r0 *
            r0;


        const double f1 =
            FOUR_PI *
            density[i] *
            r1 *
            r1;


        innerCharge[i] =
            innerCharge[i - 1] +
            0.5 *
            (f0 + f1) *
            (r1 - r0);
    }


    for (std::size_t i = N - 1;
         i > 0;
         --i)
    {
        const double r0 =
            result.r[i - 1];


        const double r1 =
            result.r[i];


        const double f0 =
            FOUR_PI *
            density[i - 1] *
            r0;


        const double f1 =
            FOUR_PI *
            density[i] *
            r1;


        outerContribution[i - 1] =
            outerContribution[i] +
            0.5 *
            (f0 + f1) *
            (r1 - r0);
    }


    for (std::size_t i = 0;
         i < N;
         ++i)
    {
        const double r =
            result.r[i];


        if (r <= 0.0)
        {
            potential[i] =
                0.0;

            continue;
        }


        potential[i] =
            innerCharge[i] / r +
            outerContribution[i];
    }
}


// ============================================================
// POTENCIAL DE INTERCAMBIO LDA
// ============================================================

void DFTHydrogen::calculateExchangePotential(
    const std::vector<double>& density,
    std::vector<double>& potential
) const
{
    const std::size_t N =
        result.r.size();


    potential.resize(N);


    const double coefficient =
        -std::cbrt(
            3.0 / PI
        );


    for (std::size_t i = 0;
         i < N;
         ++i)
    {
        const double rho =
            std::max(
                density[i],
                0.0
            );


        if (rho <= 0.0)
        {
            potential[i] =
                0.0;

            continue;
        }


        potential[i] =
            coefficient *
            std::cbrt(rho);
    }
}


// ============================================================
// POTENCIAL EFECTIVO KOHN-SHAM
// ============================================================

void DFTHydrogen::calculateEffectivePotential(
    const std::vector<double>& hartreePotential,
    const std::vector<double>& exchangePotential,
    std::vector<double>& effectivePotential
) const
{
    const std::size_t N =
        result.r.size();


    effectivePotential.resize(N);


    for (std::size_t i = 0;
         i < N;
         ++i)
    {
        const double r =
            result.r[i];


        const double nuclearPotential =
            -1.0 / r;


        effectivePotential[i] =
            nuclearPotential +
            hartreePotential[i] +
            exchangePotential[i];
    }
}


// ============================================================
// SOLVER KOHN-SHAM
// ============================================================

double DFTHydrogen::solveKohnShamOrbital(
    const std::vector<double>& effectivePotential,
    std::vector<double>& orbital
) const
{
    const std::size_t N =
        result.r.size();


    if (N < 3)
    {
        throw std::runtime_error(
            "DFTHydrogen: malla insuficiente "
            "para resolver Kohn-Sham."
        );
    }


    std::vector<double> current =
        orbital;


    std::vector<double> next(
        N,
        0.0
    );


    normalizeOrbital(
        current
    );


    const double dr =
        result.r[1] -
        result.r[0];


    const double dt =
        std::min(
            0.05,
            0.15 * dr * dr
        );


    double eigenvalue =
        -0.5;


    constexpr std::size_t MAX_ITERATIONS =
        4000;


    for (std::size_t iteration = 0;
         iteration < MAX_ITERATIONS;
         ++iteration)
    {
        std::vector<double> hamiltonian(
            N,
            0.0
        );


        hamiltonian[0] =
            0.0;


        hamiltonian[N - 1] =
            0.0;


        for (std::size_t i = 1;
             i < N - 1;
             ++i)
        {
            const double secondDerivative =
                (
                    current[i + 1]
                    - 2.0 * current[i]
                    + current[i - 1]
                ) /
                (dr * dr);


            hamiltonian[i] =
                -0.5 *
                secondDerivative +
                effectivePotential[i] *
                current[i];
        }


        const double denominator =
            dot(
                current,
                current
            );


        if (denominator <=
            std::numeric_limits<double>::epsilon())
        {
            throw std::runtime_error(
                "DFTHydrogen: orbital degenerado "
                "durante Kohn-Sham."
            );
        }


        eigenvalue =
            dot(
                current,
                hamiltonian
            ) /
            denominator;


        for (std::size_t i = 1;
             i < N - 1;
             ++i)
        {
            next[i] =
                current[i] -
                dt *
                (
                    hamiltonian[i] -
                    eigenvalue *
                    current[i]
                );
        }


        next[0] =
            0.0;


        next[N - 1] =
            0.0;


        normalizeOrbital(
            next
        );


        double difference =
            0.0;


        for (std::size_t i = 0;
             i < N;
             ++i)
        {
            const double d =
                next[i] -
                current[i];


            difference +=
                d * d;
        }


        current.swap(
            next
        );


        if (difference < 1.0e-12)
        {
            break;
        }
    }


    orbital =
        current;


    std::vector<double> hamiltonian(
        N,
        0.0
    );


    hamiltonian[0] =
        0.0;


    hamiltonian[N - 1] =
        0.0;


    for (std::size_t i = 1;
         i < N - 1;
         ++i)
    {
        const double secondDerivative =
            (
                orbital[i + 1]
                - 2.0 * orbital[i]
                + orbital[i - 1]
            ) /
            (dr * dr);


        hamiltonian[i] =
            -0.5 *
            secondDerivative +
            effectivePotential[i] *
            orbital[i];
    }


    const double denominator =
        dot(
            orbital,
            orbital
        );


    eigenvalue =
        dot(
            orbital,
            hamiltonian
        ) /
        denominator;


    return eigenvalue;
}


// ============================================================
// ENERGÍA CINÉTICA
// ============================================================

double DFTHydrogen::calculateKineticEnergy(
    const std::vector<double>& orbital
) const
{
    const std::size_t N =
        result.r.size();


    const double dr =
        result.r[1] -
        result.r[0];


    double integral =
        0.0;


    for (std::size_t i = 1;
         i < N;
         ++i)
    {
        const double derivative =
            (
                orbital[i] -
                orbital[i - 1]
            ) /
            dr;


        integral +=
            derivative *
            derivative *
            dr;
    }


    return 0.5 *
           integral;
}


// ============================================================
// ENERGÍA NÚCLEO-ELECTRÓN
// ============================================================

double DFTHydrogen::calculateNuclearAttractionEnergy(
    const std::vector<double>& density
) const
{
    const std::size_t N =
        result.r.size();


    std::vector<double> integrand(
        N,
        0.0
    );


    for (std::size_t i = 0;
         i < N;
         ++i)
    {
        const double r =
            result.r[i];


        integrand[i] =
            -FOUR_PI *
            density[i] *
            r;
    }


    return integrate(
        integrand
    );
}


// ============================================================
// ENERGÍA DE HARTREE
// ============================================================

double DFTHydrogen::calculateHartreeEnergy(
    const std::vector<double>& density,
    const std::vector<double>& hartreePotential
) const
{
    const std::size_t N =
        result.r.size();


    std::vector<double> integrand(
        N,
        0.0
    );


    for (std::size_t i = 0;
         i < N;
         ++i)
    {
        const double r =
            result.r[i];


        integrand[i] =
            0.5 *
            FOUR_PI *
            density[i] *
            hartreePotential[i] *
            r *
            r;
    }


    return integrate(
        integrand
    );
}


// ============================================================
// ENERGÍA DE INTERCAMBIO LDA
// ============================================================

double DFTHydrogen::calculateExchangeEnergy(
    const std::vector<double>& density
) const
{
    const std::size_t N =
        result.r.size();


    const double coefficient =
        -0.75 *
        std::cbrt(
            3.0 / PI
        );


    std::vector<double> integrand(
        N,
        0.0
    );


    for (std::size_t i = 0;
         i < N;
         ++i)
    {
        const double rho =
            std::max(
                density[i],
                0.0
            );


        const double r =
            result.r[i];


        integrand[i] =
            coefficient *
            std::pow(
                rho,
                4.0 / 3.0
            ) *
            FOUR_PI *
            r *
            r;
    }


    return integrate(
        integrand
    );
}


// ============================================================
// ENERGÍA TOTAL
// ============================================================

double DFTHydrogen::calculateTotalEnergy(
    double kineticEnergy,
    double nuclearAttractionEnergy,
    double hartreeEnergy,
    double exchangeEnergy
) const
{
    return
        kineticEnergy +
        nuclearAttractionEnergy +
        hartreeEnergy +
        exchangeEnergy;
}


// ============================================================
// DIFERENCIA DE DENSIDAD
// ============================================================

double DFTHydrogen::calculateDensityDifference(
    const std::vector<double>& oldDensity,
    const std::vector<double>& newDensity
) const
{
    if (oldDensity.size() !=
        newDensity.size())
    {
        throw std::invalid_argument(
            "DFTHydrogen: densidades de tamanos diferentes."
        );
    }


    if (oldDensity.empty())
    {
        return 0.0;
    }


    double maximumDifference =
        0.0;


    for (std::size_t i = 0;
         i < oldDensity.size();
         ++i)
    {
        maximumDifference =
            std::max(
                maximumDifference,
                std::abs(
                    oldDensity[i] -
                    newDensity[i]
                )
            );
    }


    return maximumDifference;
}


// ============================================================
// MEZCLA DE DENSIDADES
// ============================================================

void DFTHydrogen::mixDensity(
    const std::vector<double>& oldDensity,
    const std::vector<double>& newDensity,
    std::vector<double>& mixedDensity
) const
{
    if (oldDensity.size() !=
        newDensity.size())
    {
        throw std::invalid_argument(
            "DFTHydrogen: densidades de tamanos diferentes."
        );
    }


    mixedDensity.resize(
        oldDensity.size()
    );


    const double alpha =
        parameters.mixing;


    for (std::size_t i = 0;
         i < oldDensity.size();
         ++i)
    {
        mixedDensity[i] =
            (1.0 - alpha) *
            oldDensity[i] +
            alpha *
            newDensity[i];
    }
}


// ============================================================
// INTERPOLACIÓN LINEAL
// ============================================================

double DFTHydrogen::interpolate(
    double x,
    const std::vector<double>& xValues,
    const std::vector<double>& yValues
) const
{
    if (xValues.empty() ||
        yValues.empty() ||
        xValues.size() != yValues.size())
    {
        throw std::invalid_argument(
            "DFTHydrogen: datos invalidos "
            "para interpolacion."
        );
    }


    if (x <= xValues.front())
    {
        return yValues.front();
    }


    if (x >= xValues.back())
    {
        return yValues.back();
    }


    const auto upper =
        std::upper_bound(
            xValues.begin(),
            xValues.end(),
            x
        );


    const std::size_t i =
        static_cast<std::size_t>(
            std::distance(
                xValues.begin(),
                upper
            )
        );


    const double x0 =
        xValues[i - 1];


    const double x1 =
        xValues[i];


    const double y0 =
        yValues[i - 1];


    const double y1 =
        yValues[i];


    const double t =
        (x - x0) /
        (x1 - x0);


    return
        y0 +
        t *
        (y1 - y0);
}


// ============================================================
// INTEGRACIÓN TRAPEZOIDAL
// ============================================================

double DFTHydrogen::integrate(
    const std::vector<double>& values
) const
{
    if (values.size() !=
        result.r.size())
    {
        throw std::invalid_argument(
            "DFTHydrogen: tamanos incompatibles "
            "en integracion."
        );
    }


    if (values.size() < 2)
    {
        return 0.0;
    }


    double integral =
        0.0;


    for (std::size_t i = 1;
         i < values.size();
         ++i)
    {
        const double dx =
            result.r[i] -
            result.r[i - 1];


        integral +=
            0.5 *
            (
                values[i - 1] +
                values[i]
            ) *
            dx;
    }


    return integral;
}


// ============================================================
// PRODUCTO INTERNO
// ============================================================

double DFTHydrogen::dot(
    const std::vector<double>& a,
    const std::vector<double>& b
) const
{
    if (a.size() !=
        b.size())
    {
        throw std::invalid_argument(
            "DFTHydrogen: vectores de tamanos diferentes."
        );
    }


    if (a.empty())
    {
        return 0.0;
    }


    double sum =
        0.0;


    for (std::size_t i = 0;
         i < a.size();
         ++i)
    {
        sum +=
            a[i] *
            b[i];
    }


    if (result.r.size() >= 2)
    {
        const double dr =
            result.r[1] -
            result.r[0];


        return sum *
               dr;
    }


    return sum;
}


// ============================================================
// NORMA
// ============================================================

double DFTHydrogen::norm(
    const std::vector<double>& values
) const
{
    const double value =
        dot(
            values,
            values
        );


    if (value <= 0.0)
    {
        return 0.0;
    }


    return std::sqrt(
        value
    );
}