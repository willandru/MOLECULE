#include "NavigationViewController.h"

#include "AtomicDFT.h"
#include "AtomicOrbital3D.h"
#include "AtomicOrbitalAngular.h"
#include "AtomicOrbitalIsosurface.h"
#include "DFTConstants.h"
#include "ElectronicConfiguration.h"
#include "PBE96.h"
#include "RadialGrid.h"

#include "OrbitalRenderer.h"

#include "Grid.h"
#include "GridRenderer.h"

#include <glad/glad.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <cstddef>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{

constexpr int WINDOW_WIDTH = 1280;
constexpr int WINDOW_HEIGHT = 720;

// Cantidad de átomos que se desean visualizar.
// 0 = todos los átomos disponibles.
// 1 = H
// 2 = H, He
// 15 = H ... P
// 30 = H ... Zn
constexpr int DISPLAY_ATOM_COUNT = 3;

// Mostrar u ocultar el grid.
// true  = grid visible
// false = grid oculto
constexpr bool SHOW_GRID = false;

constexpr float PERIODIC_TABLE_SPACING_X = 5.0f;
constexpr float PERIODIC_TABLE_SPACING_Y = 7.0f;

constexpr std::size_t ORBITAL_GRID_POINTS = 61;
constexpr double ORBITAL_EXTENT = 8.0;
constexpr double ORBITAL_ISOVALUE = 0.02;

constexpr float ORBITAL_SCALE = 0.45f;

struct AtomPlacement
{
    int atomicNumber;
    int period;
    int group;
};

const std::vector<AtomPlacement> ATOM_PLACEMENTS =
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

glm::vec3 calculateAtomPosition(const AtomPlacement& placement)
{
    constexpr float CENTER_GROUP = 9.5f;
    constexpr float CENTER_PERIOD = 2.5f;

    const float x =
        (static_cast<float>(placement.group) - CENTER_GROUP) *
        PERIODIC_TABLE_SPACING_X;

    const float y =
        (CENTER_PERIOD - static_cast<float>(placement.period)) *
        PERIODIC_TABLE_SPACING_Y;

    return glm::vec3(x, y, 0.0f);
}

struct VisualOrbital
{
    AtomicOrbitalIsosurface::Result surface;
    glm::vec3 position = glm::vec3(0.0f);
    std::string atomSymbol;
    std::string orbitalName;
};

std::vector<VisualOrbital> buildAtomOrbitals(
    int Z,
    const RadialGrid& radialGrid,
    const XCFunctional& functional,
    const glm::vec3& atomPosition)
{
    const AtomicConfiguration configuration =
        getAtomicConfiguration(Z);

    const AtomicResult atom =
        solveAtom(radialGrid, configuration, functional);

    if (!atom.scf.converged)
    {
        throw std::runtime_error(
            "SCF no convergio para " + atom.symbol
        );
    }

    std::cout
        << "Atom: " << atom.symbol
        << " | Z = " << atom.Z
        << " | Electrons = " << atom.electrons
        << " | Converged = "
        << (atom.scf.converged ? "Yes" : "No")
        << '\n';

    AtomicOrbital3D orbital3D;
    AtomicOrbitalIsosurface isosurface;

    std::vector<VisualOrbital> result;

    for (const AtomicOrbital& orbital : atom.scf.orbitals)
    {
        const int angularCount =
            AtomicOrbitalAngular::orbitalCount(orbital.l);

        for (int angularIndex = 0;
             angularIndex < angularCount;
             ++angularIndex)
        {
            const AtomicOrbitalAngularType angularType =
                AtomicOrbitalAngular::typeFromQuantumNumbers(
                    orbital.l,
                    angularIndex
                );

            const std::string angularName =
                AtomicOrbitalAngular::name(angularType);

            const AtomicOrbital3D::Grid grid =
                orbital3D.sample(
                    radialGrid,
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

            visual.surface = std::move(surface);
            visual.position = atomPosition;
            visual.atomSymbol = atom.symbol;

            visual.orbitalName =
                std::to_string(orbital.n) + angularName;

            std::cout
                << "  "
                << visual.atomSymbol
                << " "
                << visual.orbitalName
                << " | + vertices = "
                << visual.surface.positive.vertices.size()
                << " | + triangles = "
                << visual.surface.positive.indices.size() / 3
                << " | - vertices = "
                << visual.surface.negative.vertices.size()
                << " | - triangles = "
                << visual.surface.negative.indices.size() / 3
                << '\n';

            result.push_back(std::move(visual));
        }
    }

    return result;
}

} // namespace

int main()
{
    try
    {
        NavigationViewController navigation(
            WINDOW_WIDTH,
            WINDOW_HEIGHT,
            "Molecule"
        );

        glEnable(GL_DEPTH_TEST);

        Grid grid;
        GridRenderer gridRenderer;
        gridRenderer.initialize(grid);

        PBE96 pbe96;

        RadialGrid radialGrid(
            DFTConstants::GRID_POINTS,
            DFTConstants::RMAX
        );

        std::vector<VisualOrbital> visualOrbitals;
        std::vector<std::unique_ptr<OrbitalRenderer>> orbitalRenderers;

        const int atomCount =
            DISPLAY_ATOM_COUNT == 0
                ? static_cast<int>(ATOM_PLACEMENTS.size())
                : DISPLAY_ATOM_COUNT;

        if (atomCount < 1 ||
            atomCount > static_cast<int>(ATOM_PLACEMENTS.size()))
        {
            throw std::runtime_error(
                "DISPLAY_ATOM_COUNT debe estar entre 1 y 30, "
                "o ser 0 para mostrar todos."
            );
        }

        for (int i = 0; i < atomCount; ++i)
        {
            const AtomPlacement& placement =
                ATOM_PLACEMENTS[static_cast<std::size_t>(i)];

            const glm::vec3 atomPosition =
                calculateAtomPosition(placement);

            std::vector<VisualOrbital> atomOrbitals =
                buildAtomOrbitals(
                    placement.atomicNumber,
                    radialGrid,
                    pbe96,
                    atomPosition
                );

            for (VisualOrbital& orbital : atomOrbitals)
            {
                visualOrbitals.push_back(
                    std::move(orbital)
                );
            }
        }

        orbitalRenderers.reserve(
            visualOrbitals.size()
        );

        for (const VisualOrbital& orbital : visualOrbitals)
        {
            auto renderer =
                std::make_unique<OrbitalRenderer>();

            renderer->initialize();

            renderer->setSurface(
                orbital.surface
            );

            renderer->setPositiveColor(
                glm::vec3(0.10f, 0.35f, 1.00f)
            );

            renderer->setNegativeColor(
                glm::vec3(1.00f, 0.15f, 0.15f)
            );

            orbitalRenderers.push_back(
                std::move(renderer)
            );
        }

        std::cout
            << "\nAtoms visualized: "
            << atomCount
            << '\n';

        std::cout
            << "Total visual orbitals: "
            << visualOrbitals.size()
            << '\n';

        while (!navigation.shouldClose())
        {
            navigation.update();

            glClear(
                GL_COLOR_BUFFER_BIT |
                GL_DEPTH_BUFFER_BIT
            );

            const glm::mat4& view =
                navigation.getViewMatrix();

            const glm::mat4& projection =
                navigation.getProjectionMatrix();

            if (SHOW_GRID)
            {
                gridRenderer.render(
                    view,
                    projection
                );
            }

            for (std::size_t i = 0;
                 i < orbitalRenderers.size();
                 ++i)
            {
                const VisualOrbital& orbital =
                    visualOrbitals[i];

                const glm::mat4 model =
                    glm::translate(
                        glm::mat4(1.0f),
                        orbital.position
                    ) *
                    glm::scale(
                        glm::mat4(1.0f),
                        glm::vec3(ORBITAL_SCALE)
                    );

                orbitalRenderers[i]->setModelMatrix(
                    model
                );

                orbitalRenderers[i]->setViewMatrix(
                    view
                );

                orbitalRenderers[i]->setProjectionMatrix(
                    projection
                );

                orbitalRenderers[i]->render();
            }

            navigation.present();
        }

        return 0;
    }
    catch (const std::exception& exception)
    {
        std::cerr
            << "\nERROR:\n"
            << exception.what()
            << '\n';

        return -1;
    }
}