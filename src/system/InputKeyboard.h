#pragma once

class Camera;


class InputKeyboard
{
public:

    explicit InputKeyboard(
        void* window
    );


    // ========================================================
    // UPDATE
    // ========================================================

    void update(
        Camera& camera,
        float deltaTime
    );


    // ========================================================
    // KEY STATE
    // ========================================================

    bool isKeyPressed(
        int key
    ) const;


    bool isKeyReleased(
        int key
    ) const;


    bool shouldClose() const;


private:

    void* window;
};