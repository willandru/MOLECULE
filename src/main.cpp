#include "system/NavigationViewController.h"

#include "geometry/grid/Grid.h"
#include "graphics/GridRenderer.h"

#include <glad/glad.h>


int main()
{
    // ========================================================
    // NAVIGATION / WINDOW
    // ========================================================

    NavigationViewController navigation(
        1280,
        720,
        "Molecule"
    );


    // ========================================================
    // OPENGL
    // ========================================================

    glEnable(GL_DEPTH_TEST);


    // ========================================================
    // GRID
    // ========================================================

    Grid grid;

    GridRenderer gridRenderer;

    gridRenderer.initialize(
        grid
    );


    // ========================================================
    // MAIN LOOP
    // ========================================================

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


        navigation.present();
    }


    return 0;
}
