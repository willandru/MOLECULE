#include "Window.h"
#include "Camera.h"
#include "InputKeyboard.h"
#include "InputMouse.h"
#include "Timer1.h"

#include "Grid.h"
#include "GridRenderer.h"

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>


int main()
{
    // ========================================================
    // SYSTEM
    // ========================================================

    Window window(
        1280,
        720,
        "MOLECULE"
    );


    Camera camera(
        glm::vec3(0.0f, 0.0f, 0.0f),
        10.0f
    );


    camera.setAspectRatio(
        window.getAspectRatio()
    );


    InputKeyboard keyboard(
        window.getHandle()
    );


    InputMouse mouse(
        window.getHandle()
    );


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

    glEnable(GL_DEPTH_TEST);


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
        // INPUT / EVENTS
        // ----------------------------------------------------

        window.pollEvents();

        mouse.update();


        // ----------------------------------------------------
        // EXIT
        // ----------------------------------------------------

        if (keyboard.shouldClose())
        {
            break;
        }


        // ----------------------------------------------------
        // CAMERA
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