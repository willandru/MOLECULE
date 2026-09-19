#pragma once

#include <glm/glm.hpp>

class Camera
{
public:
    Camera(
        const glm::vec3& target = glm::vec3(0.0f),
        float distance = 10.0f,
        float fov = 45.0f,
        float nearPlane = 0.01f,
        float farPlane = 1000.0f
    );

    void moveForward(float amount);
    void moveRight(float amount);
    void moveVertical(float amount);

    void orbit(float deltaX, float deltaY);
    void setOrientation(float yaw, float pitch);

    void zoom(float delta);

    void setTarget(const glm::vec3& target);
    void setDistance(float distance);

    void setFOV(float fov);
    void setNearPlane(float nearPlane);
    void setFarPlane(float farPlane);
    void setAspectRatio(float aspectRatio);

    glm::mat4 getViewMatrix() const;
    glm::mat4 getProjectionMatrix() const;

    const glm::vec3& getPosition() const;
    const glm::vec3& getTarget() const;

    float getDistance() const;
    float getFOV() const;
    float getNearPlane() const;
    float getFarPlane() const;
    float getAspectRatio() const;

private:
    void updatePosition();

private:
    glm::vec3 position;
    glm::vec3 target;

    float distance;

    float yaw;
    float pitch;

    float fov;
    float nearPlane;
    float farPlane;

    float aspectRatio;
};