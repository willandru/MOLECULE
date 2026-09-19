#pragma once

#include <glm/glm.hpp>

class Camera;
class InputKeyboard;
class InputMouse;
class Timer1;
class Window;

class NavigationViewController
{
public:
    NavigationViewController(
        int width,
        int height,
        const char* title
    );

    ~NavigationViewController();

    NavigationViewController(
        const NavigationViewController&
    ) = delete;

    NavigationViewController& operator=(
        const NavigationViewController&
    ) = delete;

    bool shouldClose() const;

    void update();

    void present();

    const glm::mat4& getViewMatrix() const;
    const glm::mat4& getProjectionMatrix() const;

private:
    Window* window;

    Camera* camera;
    InputKeyboard* keyboard;
    InputMouse* mouse;
    Timer1* timer;

    glm::mat4 viewMatrix;
    glm::mat4 projectionMatrix;
};