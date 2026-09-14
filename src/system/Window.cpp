#include "Window.h"

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <stdexcept>


Window::Window(
    int width,
    int height,
    const std::string& title,
    bool fullscreen
)
    : window(nullptr),
      width(width),
      height(height),
      title(title)
{
    // ========================================================
    // GLFW
    // ========================================================

    if (!glfwInit())
    {
        throw std::runtime_error(
            "Failed to initialize GLFW."
        );
    }


    // ========================================================
    // OPENGL
    // ========================================================

    glfwWindowHint(
        GLFW_CONTEXT_VERSION_MAJOR,
        3
    );

    glfwWindowHint(
        GLFW_CONTEXT_VERSION_MINOR,
        3
    );

    glfwWindowHint(
        GLFW_OPENGL_PROFILE,
        GLFW_OPENGL_CORE_PROFILE
    );

#ifdef __APPLE__

    glfwWindowHint(
        GLFW_OPENGL_FORWARD_COMPAT,
        GL_TRUE
    );

#endif


    // ========================================================
    // WINDOW
    // ========================================================

    GLFWmonitor* monitor = nullptr;

    if (fullscreen)
    {
        monitor =
            glfwGetPrimaryMonitor();

        if (monitor == nullptr)
        {
            glfwTerminate();

            throw std::runtime_error(
                "Failed to find primary monitor."
            );
        }


        const GLFWvidmode* videoMode =
            glfwGetVideoMode(monitor);

        if (videoMode == nullptr)
        {
            glfwTerminate();

            throw std::runtime_error(
                "Failed to obtain monitor video mode."
            );
        }


        this->width =
            videoMode->width;

        this->height =
            videoMode->height;


        window =
            glfwCreateWindow(
                this->width,
                this->height,
                title.c_str(),
                monitor,
                nullptr
            );
    }
    else
    {
        window =
            glfwCreateWindow(
                width,
                height,
                title.c_str(),
                nullptr,
                nullptr
            );
    }


    // ========================================================
    // WINDOW CREATION CHECK
    // ========================================================

    if (!window)
    {
        glfwTerminate();

        throw std::runtime_error(
            "Failed to create GLFW window."
        );
    }


    // ========================================================
    // OPENGL CONTEXT
    // ========================================================

    glfwMakeContextCurrent(window);


    // ========================================================
    // GLAD
    // ========================================================

    if (!gladLoadGLLoader(
        reinterpret_cast<GLADloadproc>(
            glfwGetProcAddress
        )
    ))
    {
        glfwDestroyWindow(window);

        glfwTerminate();

        throw std::runtime_error(
            "Failed to initialize GLAD."
        );
    }


    // ========================================================
    // VIEWPORT
    // ========================================================

    int framebufferWidth;
    int framebufferHeight;

    glfwGetFramebufferSize(
        window,
        &framebufferWidth,
        &framebufferHeight
    );

    glViewport(
        0,
        0,
        framebufferWidth,
        framebufferHeight
    );


    // ========================================================
    // VSYNC
    // ========================================================

    glfwSwapInterval(1);
}


// ============================================================
// DESTRUCTOR
// ============================================================

Window::~Window()
{
    if (window)
    {
        glfwDestroyWindow(window);

        window = nullptr;
    }

    glfwTerminate();
}


// ============================================================
// WINDOW STATE
// ============================================================

bool Window::shouldClose() const
{
    return glfwWindowShouldClose(
        window
    );
}


// ============================================================
// BUFFER / EVENTS
// ============================================================

void Window::swapBuffers()
{
    glfwSwapBuffers(
        window
    );
}


void Window::pollEvents()
{
    glfwPollEvents();
}


// ============================================================
// SIZE
// ============================================================

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


// ============================================================
// HANDLE
// ============================================================

GLFWwindow* Window::getHandle() const
{
    return window;
}