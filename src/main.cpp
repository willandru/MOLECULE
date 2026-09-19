#include "system/NavigationViewController.h"

#include "dft/AtomicDFT.h"
#include "dft/DFTConstants.h"
#include "dft/ElectronicConfiguration.h"
#include "dft/PZ81.h"
#include "dft/RadialGrid.h"

#include "graphics/AtomRenderer.h"
#include "geometry/grid/Grid.h"
#include "graphics/GridRenderer.h"

#include <glad/glad.h>
#include <glm/glm.hpp>

#include <iostream>

int main()
{
    NavigationViewController navigation(
        1280,
        720,
        "Molecule"
    );

    glEnable(GL_DEPTH_TEST);

    Grid grid;

    GridRenderer gridRenderer;

    gridRenderer.initialize(
        grid
    );

    PZ81 pz81;

    RadialGrid radialGrid(
        DFTConstants::GRID_POINTS,
        DFTConstants::RMAX
    );

    const AtomicConfiguration configuration =
        getAtomicConfiguration(1);

    const AtomicResult hydrogen =
        solveAtom(
            radialGrid,
            configuration,
            pz81
        );

    std::cout
        << "Atom: "
        << hydrogen.symbol
        << " | Z = "
        << hydrogen.Z
        << " | Electrons = "
        << hydrogen.electrons
        << " | Converged = "
        << (hydrogen.scf.converged ? "Yes" : "No")
        << '\n';

    AtomRenderer atomRenderer;

    atomRenderer.initialize(
        radialGrid,
        hydrogen,
        glm::vec3(0.0f, 0.0f, 0.0f),
        0.02f
    );

    while (!navigation.shouldClose())
    {
        navigation.update();

        glClear(
            GL_COLOR_BUFFER_BIT |
            GL_DEPTH_BUFFER_BIT
        );

        gridRenderer.render(
            navigation.getViewMatrix(),
            navigation.getProjectionMatrix()
        );

        atomRenderer.render(
            navigation.getViewMatrix(),
            navigation.getProjectionMatrix()
        );

        navigation.present();
    }

    return 0;
}