#include "Window.h"
#include "Camera.h"
#include "InputKeyboard.h"
#include "InputMouse.h"
#include "Timer1.h"

#include "Grid.h"
#include "GridRenderer.h"

#include "HydrogenAtom.h"
#include "HydrogenRenderer.h"

#include <glad/glad.h>

#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>


int main()
{
    Window window(
        1280,
        720,
        "HYDROGEN DFT",
        true
    );


    // ========================================================
    // CAMERA
    // ========================================================

    Camera camera(
        glm::vec3(0.0f),
        70.0f
    );

    camera.setAspectRatio(
        window.getAspectRatio()
    );

    camera.setOrientation(
        0.0f,
        glm::radians(-70.0f)
    );


    // ========================================================
    // INPUT
    // ========================================================

    InputKeyboard keyboard(
        window.getHandle()
    );

    InputMouse mouse(
        window.getHandle()
    );


    // ========================================================
    // TIMER
    // ========================================================

    Timer1 timer;


    // ========================================================
    // GRID
    // ========================================================

    Grid grid;

    GridRenderer gridRenderer;

    gridRenderer.initialize(
        grid
    );


    // ========================================================
    // HYDROGEN
    // ========================================================

    HydrogenAtom hydrogen(
        glm::vec3(0.0f)
    );

    hydrogen.calculateDFT();


    // ========================================================
    // HYDROGEN RENDERER
    // ========================================================

    HydrogenRenderer hydrogenRenderer;

    if (!hydrogenRenderer.initialize())
    {
        return -1;
    }

    hydrogenRenderer.setPosition(
        hydrogen.getPosition()
    );

    hydrogenRenderer.setScale(
        1.0f
    );

    hydrogenRenderer.setIsovalue(
        0.02
    );

    hydrogenRenderer.buildGeometry(
        hydrogen
    );


    // ========================================================
    // DFT INFORMATION
    // ========================================================

    const auto& dft =
        hydrogen.getDFTResult();

    if (!hydrogen.isDFTConverged())
    {
        return -1;
    }


    // ========================================================
    // OPENGL
    // ========================================================

    glEnable(GL_DEPTH_TEST);


    // ========================================================
    // MAIN LOOP
    // ========================================================

    while (!window.shouldClose())
    {
        timer.update();

        window.pollEvents();


        // ----------------------------------------------------
        // CLOSE
        // ----------------------------------------------------

        if (keyboard.shouldClose())
        {
            break;
        }


        // ----------------------------------------------------
        // KEYBOARD
        // ----------------------------------------------------

        keyboard.update(
            camera,
            timer.getDeltaTime()
        );


        // ----------------------------------------------------
        // MOUSE
        // ----------------------------------------------------

        mouse.update();


        if (mouse.isMiddleButtonPressed())
        {
            const glm::vec2& delta =
                mouse.getDelta();

            camera.orbit(
                delta.x,
                delta.y
            );
        }


        // ----------------------------------------------------
        // ZOOM
        // ----------------------------------------------------

        const float scroll =
            mouse.getScrollDelta();


        if (scroll != 0.0f)
        {
            camera.zoom(
                scroll
            );
        }


        mouse.clearScrollDelta();


        // ====================================================
        // CLEAR
        // ====================================================

        glClear(
            GL_COLOR_BUFFER_BIT |
            GL_DEPTH_BUFFER_BIT
        );


        // ====================================================
        // GRID
        // ====================================================

        gridRenderer.render(
            camera.getViewMatrix(),
            camera.getProjectionMatrix()
        );


        // ====================================================
        // HYDROGEN
        // ====================================================

        hydrogenRenderer.render(
            camera.getViewMatrix(),
            camera.getProjectionMatrix()
        );


        // ====================================================
        // PRESENT
        // ====================================================

        window.swapBuffers();
    }


    return 0;
}