#include "AtomBuilder.h"

#include "AtomicDFT.h"
#include "AtomicOrbital3D.h"
#include "AtomicOrbitalAngular.h"
#include "ElectronicConfiguration.h"

#include <glm/gtc/matrix_transform.hpp>

#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

const std::vector<AtomBuilder::AtomPlacement>
    AtomBuilder::ATOM_PLACEMENTS =
{
    // Periodo 1
    {1, 1, 1},
    {2, 1, 18},

    // Periodo 2
    {3, 2, 1},
    {4, 2, 2},
    {5, 2, 13},
    {6, 2, 14},
    {7, 2, 15},
    {8, 2, 16},
    {9, 2, 17},
    {10, 2, 18},

    // Periodo 3
    {11, 3, 1},
    {12, 3, 2},
    {13, 3, 13},
    {14, 3, 14},
    {15, 3, 15},
    {16, 3, 16},
    {17, 3, 17},
    {18, 3, 18},

    // Periodo 4
    {19, 4, 1},
    {20, 4, 2},
    {21, 4, 3},
    {22, 4, 4},
    {23, 4, 5},
    {24, 4, 6},
    {25, 4, 7},
    {26, 4, 8},
    {27, 4, 9},
    {28, 4, 10},
    {29, 4, 11},
    {30, 4, 12}
};

AtomBuilder::AtomBuilder(
    const RadialGrid& radialGrid,
    const XCFunctional& functional,
    int atomCount,
    float spacingX,
    float spacingY
)
    : radialGrid_(radialGrid),
      functional_(functional),
      atomCount_(atomCount),
      spacingX_(spacingX),
      spacingY_(spacingY)
{
    const int availableAtomCount =
        static_cast<int>(
            ATOM_PLACEMENTS.size()
        );

    if (atomCount_ < 0 ||
        atomCount_ > availableAtomCount)
    {
        throw std::runtime_error(
            "AtomBuilder: atomCount debe estar entre 0 y " +
            std::to_string(availableAtomCount) +
            "."
        );
    }

    if (spacingX_ <= 0.0f)
    {
        throw std::runtime_error(
            "AtomBuilder: spacingX debe ser mayor que cero."
        );
    }

    if (spacingY_ <= 0.0f)
    {
        throw std::runtime_error(
            "AtomBuilder: spacingY debe ser mayor que cero."
        );
    }

    const int count =
        atomCount_ == 0
            ? availableAtomCount
            : atomCount_;

    visualOrbitals_.reserve(
        static_cast<std::size_t>(count)
    );

    for (int i = 0; i < count; ++i)
    {
        const AtomPlacement& placement =
            ATOM_PLACEMENTS[
                static_cast<std::size_t>(i)
            ];

        const glm::vec3 atomPosition =
            calculateAtomPosition(
                placement
            );

        std::vector<VisualOrbital>
            atomOrbitals =
                buildAtomOrbitals(
                    placement.atomicNumber,
                    atomPosition
                );

        for (VisualOrbital& orbital :
             atomOrbitals)
        {
            visualOrbitals_.push_back(
                std::move(orbital)
            );
        }
    }
}

glm::vec3 AtomBuilder::calculateAtomPosition(
    const AtomPlacement& placement
) const
{
    constexpr float CENTER_GROUP = 9.5f;
    constexpr float CENTER_PERIOD = 2.5f;

    const float x =
        (
            static_cast<float>(
                placement.group
            ) -
            CENTER_GROUP
        ) *
        spacingX_;

    const float y =
        (
            CENTER_PERIOD -
            static_cast<float>(
                placement.period
            )
        ) *
        spacingY_;

    return glm::vec3(
        x,
        y,
        0.0f
    );
}

void AtomBuilder::initializeRenderer(
    VisualOrbital& orbital
)
{
    auto renderer =
        std::make_unique<OrbitalRenderer>();

    renderer->initialize();

    renderer->setSurface(
        orbital.surface
    );

    renderer->setPositiveColor(
        glm::vec3(
            0.10f,
            0.35f,
            1.00f
        )
    );

    renderer->setNegativeColor(
        glm::vec3(
            1.00f,
            0.15f,
            0.15f
        )
    );

    orbital.renderer =
        std::move(renderer);
}

std::vector<AtomBuilder::VisualOrbital>
AtomBuilder::buildAtomOrbitals(
    int atomicNumber,
    const glm::vec3& atomPosition
) const
{
    const AtomicConfiguration configuration =
        getAtomicConfiguration(
            atomicNumber
        );

    const AtomicResult atom =
        solveAtom(
            radialGrid_,
            configuration,
            functional_
        );

    if (!atom.scf.converged)
    {
        throw std::runtime_error(
            "SCF no convergio para " +
            atom.symbol
        );
    }

    std::cout
        << "Atom: "
        << atom.symbol
        << " | Z = "
        << atom.Z
        << " | Electrons = "
        << atom.electrons
        << " | Converged = "
        << (
            atom.scf.converged
                ? "Yes"
                : "No"
        )
        << '\n';

    AtomicOrbital3D orbital3D;
    AtomicOrbitalIsosurface isosurface;

    std::vector<VisualOrbital> result;

    for (const AtomicOrbital& orbital :
         atom.scf.orbitals)
    {
        const int angularCount =
            AtomicOrbitalAngular::orbitalCount(
                orbital.l
            );

        for (int angularIndex = 0;
             angularIndex < angularCount;
             ++angularIndex)
        {
            const AtomicOrbitalAngularType
                angularType =
                    AtomicOrbitalAngular::
                        typeFromQuantumNumbers(
                            orbital.l,
                            angularIndex
                        );

            const std::string angularName =
                AtomicOrbitalAngular::name(
                    angularType
                );

            const AtomicOrbital3D::Grid grid =
                orbital3D.sample(
                    radialGrid_,
                    orbital,
                    angularType,
                    ORBITAL_GRID_POINTS,
                    ORBITAL_EXTENT
                );

            AtomicOrbitalIsosurface::Result surface =
                isosurface.generate(
                    grid,
                    ORBITAL_ISOVALUE
                );

            if (surface.positive.vertices.empty() &&
                surface.negative.vertices.empty())
            {
                continue;
            }

            VisualOrbital visual;

            visual.surface =
                std::move(surface);

            visual.position =
                atomPosition;

            visual.atomSymbol =
                atom.symbol;

            visual.orbitalName =
                std::to_string(
                    orbital.n
                ) +
                angularName;

            initializeRenderer(
                visual
            );

            std::cout
                << "  "
                << visual.atomSymbol
                << " "
                << visual.orbitalName
                << " | + vertices = "
                << visual.surface
                       .positive
                       .vertices
                       .size()
                << " | + triangles = "
                << visual.surface
                       .positive
                       .indices
                       .size() / 3
                << " | - vertices = "
                << visual.surface
                       .negative
                       .vertices
                       .size()
                << " | - triangles = "
                << visual.surface
                       .negative
                       .indices
                       .size() / 3
                << '\n';

            result.push_back(
                std::move(visual)
            );
        }
    }

    return result;
}

std::size_t AtomBuilder::orbitalCount() const
{
    return visualOrbitals_.size();
}

void AtomBuilder::render(
    const glm::mat4& view,
    const glm::mat4& projection,
    float scale
)
{
    for (VisualOrbital& orbital :
         visualOrbitals_)
    {
        const glm::mat4 model =
            glm::translate(
                glm::mat4(1.0f),
                orbital.position
            ) *
            glm::scale(
                glm::mat4(1.0f),
                glm::vec3(scale)
            );

        orbital.renderer->setModelMatrix(
            model
        );

        orbital.renderer->setViewMatrix(
            view
        );

        orbital.renderer->setProjectionMatrix(
            projection
        );

        orbital.renderer->render();
    }
}