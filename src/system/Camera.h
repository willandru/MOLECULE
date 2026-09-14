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


    // ========================================================
    // CAMERA MOVEMENT
    // ========================================================

    void moveForward(
        float amount
    );

    void moveRight(
        float amount
    );


    // ========================================================
    // CAMERA ORBIT
    // ========================================================

    void orbit(
        float deltaX,
        float deltaY
    );

    void setOrientation(
        float yaw,
        float pitch
    );


    // ========================================================
    // CAMERA ZOOM
    // ========================================================

    void zoom(
        float delta
    );


    // ========================================================
    // TARGET / DISTANCE
    // ========================================================

    void setTarget(
        const glm::vec3& target
    );

    void setDistance(
        float distance
    );


    // ========================================================
    // PROJECTION
    // ========================================================

    void setFOV(
        float fov
    );

    void setNearPlane(
        float nearPlane
    );

    void setFarPlane(
        float farPlane
    );

    void setAspectRatio(
        float aspectRatio
    );


    // ========================================================
    // MATRICES
    // ========================================================

    glm::mat4 getViewMatrix() const;

    glm::mat4 getProjectionMatrix() const;


    // ========================================================
    // GETTERS
    // ========================================================

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

    // ========================================================
    // TRANSFORM
    // ========================================================

    glm::vec3 position;

    glm::vec3 target;


    // ========================================================
    // ORBIT
    // ========================================================

    float distance;

    float yaw;

    float pitch;


    // ========================================================
    // PROJECTION
    // ========================================================

    float fov;

    float nearPlane;

    float farPlane;

    float aspectRatio;
};