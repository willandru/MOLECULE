#pragma once

#include <glm/glm.hpp>

class InputMouse
{
public:
    explicit InputMouse(void* window);

    void update();

    bool isMiddleButtonPressed() const;

    const glm::vec2& getPosition() const;
    const glm::vec2& getDelta() const;
    float getScrollDelta() const;

private:
    void* window;

    glm::vec2 position;
    glm::vec2 previousPosition;
    glm::vec2 delta;

    float scrollDelta;

    bool middleButtonPressed;
};