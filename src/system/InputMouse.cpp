#include "InputMouse.h"

#include <GLFW/glfw3.h>


InputMouse::InputMouse(void* window)
    : window(window),
      position(0.0f),
      previousPosition(0.0f),
      delta(0.0f),
      scrollDelta(0.0f),
      middleButtonPressed(false)
{
}


void InputMouse::update()
{
    GLFWwindow* glfwWindow =
        static_cast<GLFWwindow*>(window);

    double mouseX;
    double mouseY;

    glfwGetCursorPos(
        glfwWindow,
        &mouseX,
        &mouseY
    );

    position = glm::vec2(
        static_cast<float>(mouseX),
        static_cast<float>(mouseY)
    );

    delta = position - previousPosition;

    previousPosition = position;

    middleButtonPressed =
        glfwGetMouseButton(
            glfwWindow,
            GLFW_MOUSE_BUTTON_MIDDLE
        ) == GLFW_PRESS;

    scrollDelta = 0.0f;
}


bool InputMouse::isMiddleButtonPressed() const
{
    return middleButtonPressed;
}


const glm::vec2& InputMouse::getPosition() const
{
    return position;
}


const glm::vec2& InputMouse::getDelta() const
{
    return delta;
}


float InputMouse::getScrollDelta() const
{
    return scrollDelta;
}