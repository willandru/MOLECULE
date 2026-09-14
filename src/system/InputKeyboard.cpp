#include "InputKeyboard.h"

#include <GLFW/glfw3.h>


InputKeyboard::InputKeyboard(void* window)
    : window(window)
{
}


bool InputKeyboard::isKeyPressed(int key) const
{
    return glfwGetKey(
        static_cast<GLFWwindow*>(window),
        key
    ) == GLFW_PRESS;
}


bool InputKeyboard::isKeyReleased(int key) const
{
    return glfwGetKey(
        static_cast<GLFWwindow*>(window),
        key
    ) == GLFW_RELEASE;
}


bool InputKeyboard::shouldClose() const
{
    return isKeyPressed(GLFW_KEY_ESCAPE);
}