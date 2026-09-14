#pragma once

#include <string>

struct GLFWwindow;

class Window
{
public:
    Window(
        int width,
        int height,
        const std::string& title,
        bool fullscreen = false
    );

    ~Window();

    Window(const Window&) = delete;
    Window& operator=(const Window&) = delete;

    bool shouldClose() const;

    void swapBuffers();
    void pollEvents();

    int getWidth() const;
    int getHeight() const;

    float getAspectRatio() const;

    GLFWwindow* getHandle() const;

private:
    GLFWwindow* window;

    int width;
    int height;

    std::string title;
};