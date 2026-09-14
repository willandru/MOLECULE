#include "Atom.h"

#include "CovalentRadiusConstants.h"


// ============================================================
// CONSTRUCTOR
// ============================================================

Atom::Atom(
    int atomicNumber,
    const glm::vec3& position
)
    : m_atomicNumber(atomicNumber),
      m_position(position),
      m_covalentRadius(
          CovalentRadius::get(atomicNumber)
      )
{
}


// ============================================================
// GETTERS
// ============================================================

int Atom::getAtomicNumber() const
{
    return m_atomicNumber;
}


const glm::vec3& Atom::getPosition() const
{
    return m_position;
}


float Atom::getCovalentRadius() const
{
    return m_covalentRadius;
}


// ============================================================
// SETTERS
// ============================================================

void Atom::setPosition(const glm::vec3& position)
{
    m_position = position;
}