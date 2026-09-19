#pragma once

#include <glm/glm.hpp>

struct GLFWwindow;

class InputMouse
{
public:
    explicit InputMouse(void* window);

    void update();

    bool isMiddleButtonPressed() const;

    const glm::vec2& getPosition() const;
    const glm::vec2& getDelta() const;

    float getScrollDelta() const;
    void clearScrollDelta();

private:
    static void scrollCallback(
        GLFWwindow* window,
        double xOffset,
        double yOffset
    );

private:
    void* window;

    glm::vec2 position;
    glm::vec2 previousPosition;
    glm::vec2 delta;

    float scrollDelta;

    bool middleButtonPressed;
    bool previousMiddleButtonPressed;
};