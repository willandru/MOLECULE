#pragma once

#include <cstddef>
#include <vector>

/*
    DFT para el átomo de hidrógeno.

    Unidades internas:
        - Distancia: Bohr
        - Energía: Hartree
        - Densidad: electrones / Bohr^3

    El cálculo es esféricamente simétrico:

        H
        ↓
        Kohn-Sham radial
        ↓
        orbital 1s
        ↓
        densidad electrónica rho(r)

    La densidad radial resultante será utilizada posteriormente
    por HydrogenRenderer para construir una isosuperficie 3D.
*/

class DFTHydrogen
{
public:

    struct Parameters
    {
        // Número de puntos del dominio radial.
        std::size_t gridPoints = 800;

        // Radio máximo del dominio, en Bohr.
        double rMax = 20.0;

        // Mezcla de densidad entre iteraciones SCF.
        double mixing = 0.35;

        // Criterio de convergencia de la densidad.
        double convergence = 1.0e-8;

        // Máximo número de iteraciones SCF.
        std::size_t maxSCFIterations = 200;
    };

    struct Result
    {
        bool converged = false;

        std::size_t iterations = 0;

        // Autovalor del orbital de Kohn-Sham.
        double eigenvalue = 0.0;

        // Componentes energéticas.
        double kineticEnergy = 0.0;
        double nuclearAttractionEnergy = 0.0;
        double hartreeEnergy = 0.0;
        double exchangeEnergy = 0.0;

        // Energía total DFT.
        double totalEnergy = 0.0;

        /*
            Coordenada radial:

                r[i]  [Bohr]

            Orbital radial reducido:

                u(r)

            relacionado con:

                psi(r) = u(r) / (sqrt(4*pi) * r)

            para un orbital s.
        */
        std::vector<double> r;
        std::vector<double> orbital;

        /*
            Densidad electrónica tridimensional esféricamente
            simétrica evaluada como función de r:

                rho(r)

            en electrones / Bohr^3.
        */
        std::vector<double> density;

        /*
            Potencial efectivo de Kohn-Sham:

                V_KS(r)
        */
        std::vector<double> effectivePotential;
    };

public:

    DFTHydrogen();

    explicit DFTHydrogen(const Parameters& parameters);

    /*
        Ejecuta todo el cálculo DFT.

        Devuelve el resultado completo del cálculo.
    */
    Result solve();

    /*
        Devuelve el último resultado calculado.
    */
    const Result& getResult() const;

private:

    Parameters parameters;
    Result result;

private:

    void buildRadialGrid();

    void initializeOrbital();

    void normalizeOrbital(
        std::vector<double>& orbital
    ) const;

    void calculateDensity(
        const std::vector<double>& orbital,
        std::vector<double>& density
    ) const;

    void calculateHartreePotential(
        const std::vector<double>& density,
        std::vector<double>& potential
    ) const;

    void calculateExchangePotential(
        const std::vector<double>& density,
        std::vector<double>& potential
    ) const;

    void calculateEffectivePotential(
        const std::vector<double>& hartreePotential,
        const std::vector<double>& exchangePotential,
        std::vector<double>& effectivePotential
    ) const;

    double solveKohnShamOrbital(
        const std::vector<double>& effectivePotential,
        std::vector<double>& orbital
    ) const;

    double calculateKineticEnergy(
        const std::vector<double>& orbital
    ) const;

    double calculateNuclearAttractionEnergy(
        const std::vector<double>& density
    ) const;

    double calculateHartreeEnergy(
        const std::vector<double>& density,
        const std::vector<double>& hartreePotential
    ) const;

    double calculateExchangeEnergy(
        const std::vector<double>& density
    ) const;

    double calculateTotalEnergy(
        double kineticEnergy,
        double nuclearAttractionEnergy,
        double hartreeEnergy,
        double exchangeEnergy
    ) const;

    double calculateDensityDifference(
        const std::vector<double>& oldDensity,
        const std::vector<double>& newDensity
    ) const;

    void mixDensity(
        const std::vector<double>& oldDensity,
        const std::vector<double>& newDensity,
        std::vector<double>& mixedDensity
    ) const;

    double interpolate(
        double x,
        const std::vector<double>& xValues,
        const std::vector<double>& yValues
    ) const;

    double integrate(
        const std::vector<double>& values
    ) const;

    double dot(
        const std::vector<double>& a,
        const std::vector<double>& b
    ) const;

    double norm(
        const std::vector<double>& values
    ) const;
};