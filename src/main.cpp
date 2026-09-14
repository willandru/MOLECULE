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
        "MOLECULE",
        true
    );


    Camera camera(
        glm::vec3(0.0f),
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
        // CAMERA SPEED
        // ----------------------------------------------------

        const float cameraSpeed =
            5.0f *
            timer.getDeltaTime();


        // ----------------------------------------------------
        // WASD
        // ----------------------------------------------------

        if (keyboard.isKeyPressed(GLFW_KEY_W))
        {
            camera.moveForward(
                cameraSpeed
            );
        }


        if (keyboard.isKeyPressed(GLFW_KEY_S))
        {
            camera.moveForward(
                -cameraSpeed
            );
        }


        if (keyboard.isKeyPressed(GLFW_KEY_A))
        {
            camera.moveRight(
                -cameraSpeed
            );
        }


        if (keyboard.isKeyPressed(GLFW_KEY_D))
        {
            camera.moveRight(
                cameraSpeed
            );
        }


        // ----------------------------------------------------
        // MIDDLE MOUSE ORBIT
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
        // MOUSE WHEEL ZOOM
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