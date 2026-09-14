#include "Window.h"

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <stdexcept>


Window::Window(int width, int height, const std::string& title)
    : window(nullptr),
      width(width),
      height(height),
      title(title)
{
    if (!glfwInit())
    {
        throw std::runtime_error("Failed to initialize GLFW.");
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    window = glfwCreateWindow(
        width,
        height,
        title.c_str(),
        nullptr,
        nullptr
    );

    if (!window)
    {
        glfwTerminate();
        throw std::runtime_error("Failed to create GLFW window.");
    }

    glfwMakeContextCurrent(window);

    if (!gladLoadGLLoader(
        reinterpret_cast<GLADloadproc>(glfwGetProcAddress)))
    {
        glfwDestroyWindow(window);
        glfwTerminate();

        throw std::runtime_error("Failed to initialize GLAD.");
    }

    glfwSwapInterval(1);
}


Window::~Window()
{
    if (window)
    {
        glfwDestroyWindow(window);
        window = nullptr;
    }

    glfwTerminate();
}


bool Window::shouldClose() const
{
    return glfwWindowShouldClose(window);
}


void Window::swapBuffers()
{
    glfwSwapBuffers(window);
}


void Window::pollEvents()
{
    glfwPollEvents();
}


int Window::getWidth() const
{
    return width;
}


int Window::getHeight() const
{
    return height;
}


float Window::getAspectRatio() const
{
    if (height == 0)
    {
        return 0.0f;
    }

    return static_cast<float>(width) /
           static_cast<float>(height);
}


GLFWwindow* Window::getHandle() const
{
    return window;
}