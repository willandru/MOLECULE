#include "DFTGeometry.h"

#include <cmath>
#include <algorithm>
#include <stdexcept>


namespace
{
    constexpr double PI =
        3.1415926535897932384626433832795;

    constexpr double EPSILON =
        1.0e-12;
}


// ================================================================
// CONSTRUCTORES
// ================================================================

DFTGeometry::DFTGeometry()
{
}


DFTGeometry::DFTGeometry(
    const std::vector<int>& nuclearCharges,
    const std::vector<glm::dvec3>& nuclearPositions
)
{
    initialize(
        nuclearCharges,
        nuclearPositions
    );
}


// ================================================================
// INICIALIZACIÓN
// ================================================================

void DFTGeometry::initialize(
    const std::vector<int>& charges,
    const std::vector<glm::dvec3>& positions
)
{
    if (charges.size() != positions.size())
    {
        throw std::invalid_argument(
            "DFTGeometry: charges and positions have different sizes."
        );
    }

    for (int charge : charges)
    {
        if (charge <= 0)
        {
            throw std::invalid_argument(
                "DFTGeometry: nuclear charge must be positive."
            );
        }
    }

    nuclearCharges =
        charges;

    nuclearPositions =
        positions;

    forces.assign(
        nuclearPositions.size(),
        glm::dvec3(0.0)
    );
}


// ================================================================
// NÚCLEOS
// ================================================================

std::size_t DFTGeometry::getNuclearCount() const
{
    return nuclearPositions.size();
}


int DFTGeometry::getNuclearCharge(
    std::size_t index
) const
{
    return nuclearCharges.at(index);
}


const glm::dvec3&
DFTGeometry::getNuclearPosition(
    std::size_t index
) const
{
    return nuclearPositions.at(index);
}


glm::dvec3&
DFTGeometry::getNuclearPosition(
    std::size_t index
)
{
    return nuclearPositions.at(index);
}


const std::vector<int>&
DFTGeometry::getNuclearCharges() const
{
    return nuclearCharges;
}


const std::vector<glm::dvec3>&
DFTGeometry::getNuclearPositions() const
{
    return nuclearPositions;
}


std::vector<glm::dvec3>&
DFTGeometry::getNuclearPositions()
{
    return nuclearPositions;
}


void DFTGeometry::setNuclearPosition(
    std::size_t index,
    const glm::dvec3& position
)
{
    if (!std::isfinite(position.x) ||
        !std::isfinite(position.y) ||
        !std::isfinite(position.z))
    {
        throw std::invalid_argument(
            "DFTGeometry: position contains a non-finite value."
        );
    }

    nuclearPositions.at(index) =
        position;
}


// ================================================================
// DISTANCIAS
// ================================================================

double DFTGeometry::calculateDistance(
    std::size_t first,
    std::size_t second
) const
{
    const glm::dvec3 difference =
        nuclearPositions.at(first) -
        nuclearPositions.at(second);

    return glm::length(
        difference
    );
}


// ================================================================
// ÁNGULOS
// ================================================================

double DFTGeometry::calculateAngle(
    std::size_t first,
    std::size_t center,
    std::size_t second
) const
{
    /*
        Para:

            first ---- center ---- second

        definimos:

            a = R_first  - R_center
            b = R_second - R_center

        y:

            cos(theta) =
            (a·b)/(|a||b|)
    */

    const glm::dvec3 a =
        nuclearPositions.at(first) -
        nuclearPositions.at(center);

    const glm::dvec3 b =
        nuclearPositions.at(second) -
        nuclearPositions.at(center);

    const double lengthA =
        glm::length(a);

    const double lengthB =
        glm::length(b);

    if (lengthA <= EPSILON ||
        lengthB <= EPSILON)
    {
        throw std::runtime_error(
            "DFTGeometry: angle contains overlapping nuclei."
        );
    }

    double cosine =
        glm::dot(a, b) /
        (lengthA * lengthB);

    cosine =
        std::clamp(
            cosine,
            -1.0,
            1.0
        );

    return std::acos(
        cosine
    );
}


// ================================================================
// FUERZAS
// ================================================================

void DFTGeometry::setForces(
    const std::vector<glm::dvec3>& newForces
)
{
    if (newForces.size() !=
        nuclearPositions.size())
    {
        throw std::invalid_argument(
            "DFTGeometry: force count does not match nuclear count."
        );
    }

    for (const glm::dvec3& force : newForces)
    {
        if (!std::isfinite(force.x) ||
            !std::isfinite(force.y) ||
            !std::isfinite(force.z))
        {
            throw std::invalid_argument(
                "DFTGeometry: force contains a non-finite value."
            );
        }
    }

    forces =
        newForces;
}


const std::vector<glm::dvec3>&
DFTGeometry::getForces() const
{
    return forces;
}


glm::dvec3 DFTGeometry::getForce(
    std::size_t index
) const
{
    return forces.at(index);
}


double DFTGeometry::getForceNorm(
    std::size_t index
) const
{
    return glm::length(
        forces.at(index)
    );
}


double DFTGeometry::getMaximumForce() const
{
    double maximum =
        0.0;

    for (const glm::dvec3& force : forces)
    {
        maximum =
            std::max(
                maximum,
                glm::length(force)
            );
    }

    return maximum;
}


// ================================================================
// MOVIMIENTO
// ================================================================

void DFTGeometry::moveNuclei(
    double step
)
{
    if (step <= 0.0 ||
        !std::isfinite(step))
    {
        throw std::invalid_argument(
            "DFTGeometry: step must be positive and finite."
        );
    }

    /*
        F = -∇E

        Para minimizar E:

            R_new = R_old + αF

        porque la fuerza apunta hacia
        la dirección de descenso de energía.
    */

    for (std::size_t i = 0;
         i < nuclearPositions.size();
         ++i)
    {
        nuclearPositions[i] +=
            step *
            forces[i];
    }
}


void DFTGeometry::moveNucleus(
    std::size_t index,
    const glm::dvec3& displacement
)
{
    if (!std::isfinite(displacement.x) ||
        !std::isfinite(displacement.y) ||
        !std::isfinite(displacement.z))
    {
        throw std::invalid_argument(
            "DFTGeometry: displacement contains a non-finite value."
        );
    }

    nuclearPositions.at(index) +=
        displacement;
}


// ================================================================
// OPTIMIZACIÓN
// ================================================================

bool DFTGeometry::optimizeStep(
    const std::vector<glm::dvec3>& newForces,
    double step,
    double forceTolerance
)
{
    if (newForces.size() !=
        nuclearPositions.size())
    {
        throw std::invalid_argument(
            "DFTGeometry: force count does not match nuclear count."
        );
    }

    if (step <= 0.0 ||
        !std::isfinite(step))
    {
        throw std::invalid_argument(
            "DFTGeometry: optimization step must be positive and finite."
        );
    }

    if (forceTolerance <= 0.0 ||
        !std::isfinite(forceTolerance))
    {
        throw std::invalid_argument(
            "DFTGeometry: force tolerance must be positive and finite."
        );
    }

    setForces(
        newForces
    );

    if (isOptimized(forceTolerance))
    {
        return true;
    }

    moveNuclei(
        step
    );

    return false;
}


bool DFTGeometry::isOptimized(
    double forceTolerance
) const
{
    if (forceTolerance <= 0.0 ||
        !std::isfinite(forceTolerance))
    {
        throw std::invalid_argument(
            "DFTGeometry: force tolerance must be positive and finite."
        );
    }

    return
        getMaximumForce() <=
        forceTolerance;
}


// ================================================================
// CENTRAR GEOMETRÍA
// ================================================================

void DFTGeometry::centerAtOrigin()
{
    if (nuclearPositions.empty())
    {
        return;
    }

    glm::dvec3 center(
        0.0
    );

    for (const glm::dvec3& position :
         nuclearPositions)
    {
        center +=
            position;
    }

    center /=
        static_cast<double>(
            nuclearPositions.size()
        );

    for (glm::dvec3& position :
         nuclearPositions)
    {
        position -=
            center;
    }
}