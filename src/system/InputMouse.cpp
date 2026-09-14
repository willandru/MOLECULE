#include "InputMouse.h"

#include <GLFW/glfw3.h>


// ============================================================
// CONSTRUCTOR
// ============================================================

InputMouse::InputMouse(
    void* window
)
    : window(window),
      position(0.0f),
      previousPosition(0.0f),
      delta(0.0f),
      scrollDelta(0.0f),
      middleButtonPressed(false)
{
    GLFWwindow* glfwWindow =
        static_cast<GLFWwindow*>(window);


    glfwSetWindowUserPointer(
        glfwWindow,
        this
    );


    glfwSetScrollCallback(
        glfwWindow,
        scrollCallback
    );
}


// ============================================================
// UPDATE
// ============================================================

void InputMouse::update()
{
    GLFWwindow* glfwWindow =
        static_cast<GLFWwindow*>(window);


    // --------------------------------------------------------
    // POSITION
    // --------------------------------------------------------

    double mouseX;
    double mouseY;


    glfwGetCursorPos(
        glfwWindow,
        &mouseX,
        &mouseY
    );


    position =
        glm::vec2(
            static_cast<float>(mouseX),
            static_cast<float>(mouseY)
        );


    // --------------------------------------------------------
    // DELTA
    // --------------------------------------------------------

    delta =
        position -
        previousPosition;


    previousPosition =
        position;


    // --------------------------------------------------------
    // MIDDLE BUTTON
    // --------------------------------------------------------

    middleButtonPressed =
        glfwGetMouseButton(
            glfwWindow,
            GLFW_MOUSE_BUTTON_MIDDLE
        ) == GLFW_PRESS;
}


// ============================================================
// SCROLL CALLBACK
// ============================================================

void InputMouse::scrollCallback(
    GLFWwindow* window,
    double,
    double yOffset
)
{
    InputMouse* mouse =
        static_cast<InputMouse*>(
            glfwGetWindowUserPointer(window)
        );


    if (mouse != nullptr)
    {
        mouse->scrollDelta +=
            static_cast<float>(yOffset);
    }
}


// ============================================================
// MIDDLE BUTTON
// ============================================================

bool InputMouse::isMiddleButtonPressed() const
{
    return middleButtonPressed;
}


// ============================================================
// POSITION
// ============================================================

const glm::vec2& InputMouse::getPosition() const
{
    return position;
}


// ============================================================
// DELTA
// ============================================================

const glm::vec2& InputMouse::getDelta() const
{
    return delta;
}


// ============================================================
// SCROLL
// ============================================================

float InputMouse::getScrollDelta() const
{
    return scrollDelta;
}


// ============================================================
// CLEAR SCROLL
// ============================================================

void InputMouse::clearScrollDelta()
{
    scrollDelta = 0.0f;
}