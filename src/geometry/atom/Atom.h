#pragma once

#include <glm/vec3.hpp>


class Atom
{
public:

    // ========================================================
    // CONSTRUCTOR
    // ========================================================

    Atom(
        int atomicNumber,
        const glm::vec3& position
    );


    // ========================================================
    // GETTERS
    // ========================================================

    int getAtomicNumber() const;

    const glm::vec3& getPosition() const;

    float getCovalentRadius() const;


    // ========================================================
    // SETTERS
    // ========================================================

    void setPosition(const glm::vec3& position);


private:

    // ========================================================
    // ATOMIC DATA
    // ========================================================

    int m_atomicNumber;

    // Position in internal geometry units:
    // 1.0 unit = 1 Å
    glm::vec3 m_position;

    // Covalent radius in Å
    float m_covalentRadius;
};