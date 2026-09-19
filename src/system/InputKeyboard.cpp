#include "InputKeyboard.h"

#include "Camera.h"

#include <GLFW/glfw3.h>


// ============================================================
// CONSTRUCTOR
// ============================================================

InputKeyboard::InputKeyboard(
    void* window
)
    : window(window)
{
}


// ============================================================
// UPDATE
// ============================================================

void InputKeyboard::update(
    Camera& camera,
    float deltaTime
)
{
    constexpr float cameraSpeed =
        5.0f;

    const float movement =
        cameraSpeed *
        deltaTime;


    // --------------------------------------------------------
    // FORWARD / BACKWARD
    // --------------------------------------------------------

    if (isKeyPressed(GLFW_KEY_W))
    {
        camera.moveForward(
            movement
        );
    }

    if (isKeyPressed(GLFW_KEY_S))
    {
        camera.moveForward(
            -movement
        );
    }


    // --------------------------------------------------------
    // LEFT / RIGHT
    // --------------------------------------------------------

    if (isKeyPressed(GLFW_KEY_A))
    {
        camera.moveRight(
            -movement
        );
    }

    if (isKeyPressed(GLFW_KEY_D))
    {
        camera.moveRight(
            movement
        );
    }


    // --------------------------------------------------------
    // VERTICAL
    // --------------------------------------------------------

    if (isKeyPressed(GLFW_KEY_Q))
    {
        camera.moveVertical(
            -movement
        );
    }

    if (isKeyPressed(GLFW_KEY_E))
    {
        camera.moveVertical(
            movement
        );
    }
}


// ============================================================
// KEY PRESSED
// ============================================================

bool InputKeyboard::isKeyPressed(
    int key
) const
{
    return glfwGetKey(
        static_cast<GLFWwindow*>(window),
        key
    ) == GLFW_PRESS;
}


// ============================================================
// KEY RELEASED
// ============================================================

bool InputKeyboard::isKeyReleased(
    int key
) const
{
    return glfwGetKey(
        static_cast<GLFWwindow*>(window),
        key
    ) == GLFW_RELEASE;
}


// ============================================================
// SHOULD CLOSE
// ============================================================

bool InputKeyboard::shouldClose() const
{
    return isKeyPressed(
        GLFW_KEY_ESCAPE
    );
}