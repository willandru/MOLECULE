#pragma once

#include <vector>
#include <cstddef>

#include <glm/glm.hpp>

class DFTGeometry
{
public:

    DFTGeometry();

    DFTGeometry(
        const std::vector<int>& nuclearCharges,
        const std::vector<glm::dvec3>& nuclearPositions
    );

    void initialize(
        const std::vector<int>& nuclearCharges,
        const std::vector<glm::dvec3>& nuclearPositions
    );

    // ------------------------------------------------------------
    // Núcleos
    // ------------------------------------------------------------

    std::size_t getNuclearCount() const;

    int getNuclearCharge(
        std::size_t index
    ) const;

    const glm::dvec3& getNuclearPosition(
        std::size_t index
    ) const;

    glm::dvec3& getNuclearPosition(
        std::size_t index
    );

    const std::vector<int>& getNuclearCharges() const;

    const std::vector<glm::dvec3>& getNuclearPositions() const;

    std::vector<glm::dvec3>& getNuclearPositions();

    void setNuclearPosition(
        std::size_t index,
        const glm::dvec3& position
    );

    // ------------------------------------------------------------
    // Distancias
    // ------------------------------------------------------------

    double calculateDistance(
        std::size_t first,
        std::size_t second
    ) const;

    // ------------------------------------------------------------
    // Ángulos
    // ------------------------------------------------------------

    double calculateAngle(
        std::size_t first,
        std::size_t center,
        std::size_t second
    ) const;

    // ------------------------------------------------------------
    // Fuerzas
    // ------------------------------------------------------------

    void setForces(
        const std::vector<glm::dvec3>& forces
    );

    const std::vector<glm::dvec3>& getForces() const;

    glm::dvec3 getForce(
        std::size_t index
    ) const;

    double getMaximumForce() const;

    double getForceNorm(
        std::size_t index
    ) const;

    // ------------------------------------------------------------
    // Movimiento
    // ------------------------------------------------------------

    void moveNuclei(
        double step
    );

    void moveNucleus(
        std::size_t index,
        const glm::dvec3& displacement
    );

    // ------------------------------------------------------------
    // Optimización
    // ------------------------------------------------------------

    bool optimizeStep(
        const std::vector<glm::dvec3>& forces,
        double step,
        double forceTolerance
    );

    bool isOptimized(
        double forceTolerance
    ) const;

    // ------------------------------------------------------------
    // Utilidades
    // ------------------------------------------------------------

    void centerAtOrigin();

private:

    std::vector<int> nuclearCharges;

    std::vector<glm::dvec3> nuclearPositions;

    std::vector<glm::dvec3> forces;
};