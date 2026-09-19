#pragma once

#include <string>

enum class AtomicOrbitalAngularType
{
    S,

    Px,
    Py,
    Pz,

    Dxy,
    Dxz,
    Dyz,
    Dx2y2,
    Dz2
};

class AtomicOrbitalAngular
{
public:

    static AtomicOrbitalAngularType typeFromQuantumNumbers(
        int l,
        int orbitalIndex
    );

    static int orbitalCount(
        int l
    );

    static std::string name(
        AtomicOrbitalAngularType type
    );

    static double evaluate(
        AtomicOrbitalAngularType type,
        double x,
        double y,
        double z
    );

    static double evaluateNormalized(
        AtomicOrbitalAngularType type,
        double x,
        double y,
        double z
    );
};