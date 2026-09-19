#include "AtomicOrbitalAngular.h"

#include <cmath>
#include <stdexcept>

namespace
{

constexpr double PI =
    3.1415926535897932384626433832795;

double radius(
    double x,
    double y,
    double z
)
{
    return std::sqrt(
        x * x +
        y * y +
        z * z
    );
}

}

AtomicOrbitalAngularType
AtomicOrbitalAngular::typeFromQuantumNumbers(
    int l,
    int orbitalIndex
)
{
    if (orbitalIndex < 0) {
        throw std::invalid_argument(
            "El indice orbital no puede ser negativo."
        );
    }

    switch (l) {

    case 0:

        if (orbitalIndex == 0) {
            return AtomicOrbitalAngularType::S;
        }

        break;

    case 1:

        switch (orbitalIndex) {

        case 0:
            return AtomicOrbitalAngularType::Px;

        case 1:
            return AtomicOrbitalAngularType::Py;

        case 2:
            return AtomicOrbitalAngularType::Pz;

        default:
            break;
        }

        break;

    case 2:

        switch (orbitalIndex) {

        case 0:
            return AtomicOrbitalAngularType::Dxy;

        case 1:
            return AtomicOrbitalAngularType::Dxz;

        case 2:
            return AtomicOrbitalAngularType::Dyz;

        case 3:
            return AtomicOrbitalAngularType::Dx2y2;

        case 4:
            return AtomicOrbitalAngularType::Dz2;

        default:
            break;
        }

        break;

    default:
        break;
    }

    throw std::invalid_argument(
        "Combinacion de l e indice orbital no soportada."
    );
}

int AtomicOrbitalAngular::orbitalCount(
    int l
)
{
    if (l < 0) {
        throw std::invalid_argument(
            "El numero cuantico l no puede ser negativo."
        );
    }

    return 2 * l + 1;
}

std::string AtomicOrbitalAngular::name(
    AtomicOrbitalAngularType type
)
{
    switch (type) {

    case AtomicOrbitalAngularType::S:
        return "s";

    case AtomicOrbitalAngularType::Px:
        return "px";

    case AtomicOrbitalAngularType::Py:
        return "py";

    case AtomicOrbitalAngularType::Pz:
        return "pz";

    case AtomicOrbitalAngularType::Dxy:
        return "dxy";

    case AtomicOrbitalAngularType::Dxz:
        return "dxz";

    case AtomicOrbitalAngularType::Dyz:
        return "dyz";

    case AtomicOrbitalAngularType::Dx2y2:
        return "dx2-y2";

    case AtomicOrbitalAngularType::Dz2:
        return "dz2";
    }

    throw std::invalid_argument(
        "Tipo de orbital angular no valido."
    );
}

double AtomicOrbitalAngular::evaluate(
    AtomicOrbitalAngularType type,
    double x,
    double y,
    double z
)
{
    const double r =
        radius(x, y, z);

    if (r <= 0.0) {

        if (type ==
            AtomicOrbitalAngularType::S) {

            return 1.0;
        }

        return 0.0;
    }

    const double r2 =
        r * r;

    switch (type) {

    case AtomicOrbitalAngularType::S:

        return 1.0;

    case AtomicOrbitalAngularType::Px:

        return x / r;

    case AtomicOrbitalAngularType::Py:

        return y / r;

    case AtomicOrbitalAngularType::Pz:

        return z / r;

    case AtomicOrbitalAngularType::Dxy:

        return (x * y) / r2;

    case AtomicOrbitalAngularType::Dxz:

        return (x * z) / r2;

    case AtomicOrbitalAngularType::Dyz:

        return (y * z) / r2;

    case AtomicOrbitalAngularType::Dx2y2:

        return (x * x - y * y) / r2;

    case AtomicOrbitalAngularType::Dz2:

        return
            (2.0 * z * z -
             x * x -
             y * y) / r2;
    }

    throw std::invalid_argument(
        "Tipo de orbital angular no valido."
    );
}

double AtomicOrbitalAngular::evaluateNormalized(
    AtomicOrbitalAngularType type,
    double x,
    double y,
    double z
)
{
    const double angular =
        evaluate(
            type,
            x,
            y,
            z
        );

    switch (type) {

    case AtomicOrbitalAngularType::S:

        return
            angular *
            std::sqrt(
                1.0 /
                (4.0 * PI)
            );

    case AtomicOrbitalAngularType::Px:
    case AtomicOrbitalAngularType::Py:
    case AtomicOrbitalAngularType::Pz:

        return
            angular *
            std::sqrt(
                3.0 /
                (4.0 * PI)
            );

    case AtomicOrbitalAngularType::Dxy:
    case AtomicOrbitalAngularType::Dxz:
    case AtomicOrbitalAngularType::Dyz:

        return
            angular *
            std::sqrt(
                15.0 /
                (4.0 * PI)
            );

    case AtomicOrbitalAngularType::Dx2y2:

        return
            angular *
            std::sqrt(
                15.0 /
                (16.0 * PI)
            );

    case AtomicOrbitalAngularType::Dz2:

        return
            angular *
            std::sqrt(
                5.0 /
                (16.0 * PI)
            );
    }

    throw std::invalid_argument(
        "Tipo de orbital angular no valido."
    );
}