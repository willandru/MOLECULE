#include "Timer1.h"

#include <GLFW/glfw3.h>


Timer1::Timer1()
    : previousTime(glfwGetTime()),
      deltaTime(0.0f),
      elapsedTime(0.0f)
{
}


void Timer1::update()
{
    double currentTime = glfwGetTime();

    deltaTime = static_cast<float>(
        currentTime - previousTime
    );

    previousTime = currentTime;

    elapsedTime += deltaTime;
}


float Timer1::getDeltaTime() const
{
    return deltaTime;
}


float Timer1::getElapsedTime() const
{
    return elapsedTime;
}