#pragma once

class Timer1
{
public:
    Timer1();

    void update();

    float getDeltaTime() const;
    float getElapsedTime() const;

private:
    double previousTime;
    float deltaTime;
    float elapsedTime;
};