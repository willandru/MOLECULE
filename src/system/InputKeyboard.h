#pragma once

class InputKeyboard
{
public:
    explicit InputKeyboard(void* window);

    bool isKeyPressed(int key) const;
    bool isKeyReleased(int key) const;

    bool shouldClose() const;

private:
    void* window;
};