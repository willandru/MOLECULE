#include "NavigationViewController.h"

#include "AtomBuilder.h"
#include "DFTConstants.h"
#include "PBE96.h"
#include "RadialGrid.h"

#include "Grid.h"
#include "GridRenderer.h"

#include <glad/glad.h>

#include <glm/glm.hpp>

#include <iostream>
#include <stdexcept>

namespace
{

constexpr int WINDOW_WIDTH = 1280;
constexpr int WINDOW_HEIGHT = 720;

// 0 = todos
// 1 = H
// 2 = H, He
// 15 = H ... P
// 30 = H ... Zn
constexpr int DISPLAY_ATOM_COUNT = 1;

// true  = mostrar grid
// false = ocultar grid
constexpr bool SHOW_GRID = false;

// Separación de la tabla periódica.
constexpr float PERIODIC_TABLE_SPACING_X = 5.0f;
constexpr float PERIODIC_TABLE_SPACING_Y = 7.0f;

// Escala visual de los orbitales.
constexpr float ORBITAL_SCALE = 0.45f;

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

        gridRenderer.initialize(
            grid
        );

        PBE96 pbe96;

        RadialGrid radialGrid(
            DFTConstants::GRID_POINTS,
            DFTConstants::RMAX
        );

        AtomBuilder atomBuilder(
            radialGrid,
            pbe96,
            DISPLAY_ATOM_COUNT,
            PERIODIC_TABLE_SPACING_X,
            PERIODIC_TABLE_SPACING_Y
        );

        std::cout
            << "\nAtoms visualized: "
            << (
                DISPLAY_ATOM_COUNT == 0
                    ? 30
                    : DISPLAY_ATOM_COUNT
            )
            << '\n';

        std::cout
            << "Total visual orbitals: "
            << atomBuilder.orbitalCount()
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

            atomBuilder.render(
                view,
                projection,
                ORBITAL_SCALE
            );

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