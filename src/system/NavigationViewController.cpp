#include "NavigationViewController.h"

#include "Camera.h"
#include "InputKeyboard.h"
#include "InputMouse.h"
#include "Timer1.h"
#include "Window.h"

#include <glad/glad.h>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <stdexcept>


// ============================================================
// CONSTRUCTOR
// ============================================================

NavigationViewController::NavigationViewController(
    int width,
    int height,
    const char* title
)
    : window(nullptr),
      camera(nullptr),
      keyboard(nullptr),
      mouse(nullptr),
      timer(nullptr),
      viewMatrix(1.0f),
      projectionMatrix(1.0f)
{
    window =
        new Window(
            width,
            height,
            title,
            true
        );


    if (window->getHandle() == nullptr)
    {
        delete window;
        window = nullptr;

        throw std::runtime_error(
            "Failed to create application window."
        );
    }


    // ========================================================
    // OPENGL VIEWPORT
    // ========================================================

    glViewport(
        0,
        0,
        window->getWidth(),
        window->getHeight()
    );


    // ========================================================
    // CAMERA
    // ========================================================

    camera =
        new Camera(
            glm::vec3(0.0f),
            95.0f,
            45.0f,
            0.01f,
            1000.0f
        );


    // Vista frontal de la tabla periódica.
    //
    // Yaw = 180°  -> cámara frente al eje Z.
    // Pitch = 0°  -> sin inclinación vertical.
    //
    camera->setOrientation(
        glm::radians(180.0f),
        glm::radians(0.0f)
    );


    camera->setAspectRatio(
        window->getAspectRatio()
    );


    // ========================================================
    // INPUT
    // ========================================================

    keyboard =
        new InputKeyboard(
            static_cast<void*>(
                window->getHandle()
            )
        );


    mouse =
        new InputMouse(
            static_cast<void*>(
                window->getHandle()
            )
        );


    // ========================================================
    // TIMER
    // ========================================================

    timer =
        new Timer1();


    // ========================================================
    // INITIAL MATRICES
    // ========================================================

    viewMatrix =
        camera->getViewMatrix();


    projectionMatrix =
        camera->getProjectionMatrix();
}


// ============================================================
// DESTRUCTOR
// ============================================================

NavigationViewController::~NavigationViewController()
{
    delete timer;
    delete mouse;
    delete keyboard;
    delete camera;
    delete window;
}


// ============================================================
// SHOULD CLOSE
// ============================================================

bool NavigationViewController::shouldClose() const
{
    if (window == nullptr)
    {
        return true;
    }


    if (window->shouldClose())
    {
        return true;
    }


    if (
        keyboard != nullptr &&
        keyboard->shouldClose()
    )
    {
        return true;
    }


    return false;
}


// ============================================================
// UPDATE
// ============================================================

void NavigationViewController::update()
{
    if (
        window == nullptr ||
        camera == nullptr ||
        keyboard == nullptr ||
        mouse == nullptr ||
        timer == nullptr
    )
    {
        return;
    }


    // ========================================================
    // TIMER
    // ========================================================

    timer->update();


    const float deltaTime =
        timer->getDeltaTime();


    // ========================================================
    // KEYBOARD
    // ========================================================

    keyboard->update(
        *camera,
        deltaTime
    );


    // ========================================================
    // MOUSE
    // ========================================================

    mouse->update();


    // ========================================================
    // ORBIT
    // ========================================================

    if (mouse->isMiddleButtonPressed())
    {
        const glm::vec2 delta =
            mouse->getDelta();


        camera->orbit(
            delta.x,
            delta.y
        );
    }


    // ========================================================
    // ZOOM
    // ========================================================

    const float scrollDelta =
        mouse->getScrollDelta();


    if (scrollDelta != 0.0f)
    {
        camera->zoom(
            scrollDelta
        );


        mouse->clearScrollDelta();
    }


    // ========================================================
    // ASPECT RATIO
    // ========================================================

    camera->setAspectRatio(
        window->getAspectRatio()
    );


    // ========================================================
    // MATRICES
    // ========================================================

    viewMatrix =
        camera->getViewMatrix();


    projectionMatrix =
        camera->getProjectionMatrix();
}


// ============================================================
// PRESENT
// ============================================================

void NavigationViewController::present()
{
    if (window == nullptr)
    {
        return;
    }


    window->swapBuffers();
    window->pollEvents();
}


// ============================================================
// VIEW MATRIX
// ============================================================

const glm::mat4&
NavigationViewController::getViewMatrix() const
{
    return viewMatrix;
}


// ============================================================
// PROJECTION MATRIX
// ============================================================

const glm::mat4&
NavigationViewController::getProjectionMatrix() const
{
    return projectionMatrix;
}