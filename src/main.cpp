#include "Window.h"
#include "Camera.h"
#include "InputKeyboard.h"
#include "InputMouse.h"
#include "Timer1.h"

#include "Grid.h"
#include "GridRenderer.h"

#include <glad/glad.h>

#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>


int main()
{
    // ========================================================
    // SYSTEM
    // ========================================================

    Window window(
        1280,
        720,
        "MOLECULE",
        true
    );


    // ========================================================
    // CAMERA
    // ========================================================

    Camera camera(
        glm::vec3(0.0f),
        10.0f
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
    // TIME
    // ========================================================

    Timer1 timer;


    // ========================================================
    // GEOMETRY
    // ========================================================

    Grid grid;


    // ========================================================
    // GRAPHICS
    // ========================================================

    GridRenderer gridRenderer;

    gridRenderer.initialize(
        grid
    );


    // ========================================================
    // OPENGL
    // ========================================================

    glEnable(
        GL_DEPTH_TEST
    );


    // ========================================================
    // MAIN LOOP
    // ========================================================

    while (!window.shouldClose())
    {
        // ----------------------------------------------------
        // TIME
        // ----------------------------------------------------

        timer.update();


        // ----------------------------------------------------
        // EVENTS
        // ----------------------------------------------------

        window.pollEvents();


        // ----------------------------------------------------
        // KEYBOARD
        // ----------------------------------------------------

        if (keyboard.shouldClose())
        {
            break;
        }


        keyboard.update(
            camera,
            timer.getDeltaTime()
        );


        // ----------------------------------------------------
        // MOUSE
        // ----------------------------------------------------

        mouse.update();


        // ----------------------------------------------------
        // CAMERA ORBIT
        // ----------------------------------------------------

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
        // CAMERA ZOOM
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


        // ----------------------------------------------------
        // RENDER
        // ----------------------------------------------------

        glClear(
            GL_COLOR_BUFFER_BIT |
            GL_DEPTH_BUFFER_BIT
        );


        gridRenderer.render(
            camera.getViewMatrix(),
            camera.getProjectionMatrix()
        );


        // ----------------------------------------------------
        // PRESENT
        // ----------------------------------------------------

        window.swapBuffers();
    }


    return 0;
}