#pragma once

#include <cstddef>
#include <vector>

#include <glm/vec3.hpp>

#include "DFTHydrogen.h"


class HydrogenAtom
{
public:

    HydrogenAtom();

    explicit HydrogenAtom(
        const glm::vec3& position
    );


    // ========================================================
    // POSITION
    // ========================================================

    void setPosition(
        const glm::vec3& position
    );

    const glm::vec3& getPosition() const;


    // ========================================================
    // ATOMIC DATA
    // ========================================================

    int getAtomicNumber() const;

    int getElectronCount() const;


    // ========================================================
    // DFT
    // ========================================================

    void calculateDFT();

    bool isDFTConverged() const;

    std::size_t getSCFIterations() const;


    // ========================================================
    // ENERGY
    // ========================================================

    double getEnergy() const;

    double getKSEigenvalue() const;

    double getKineticEnergy() const;

    double getNuclearAttractionEnergy() const;

    double getHartreeEnergy() const;

    double getExchangeEnergy() const;


    // ========================================================
    // DFT DATA
    // ========================================================

    const std::vector<double>& getRadialGrid() const;

    const std::vector<double>& getOrbital() const;

    const std::vector<double>& getDensity() const;

    const std::vector<double>& getEffectivePotential() const;


    const DFTHydrogen::Result& getDFTResult() const;


private:

    static constexpr int ATOMIC_NUMBER = 1;

    static constexpr int ELECTRON_COUNT = 1;


    glm::vec3 position;

    DFTHydrogen dft;

    bool dftCalculated;
};