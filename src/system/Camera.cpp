#include "Camera.h"

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/constants.hpp>

#include <algorithm>
#include <cmath>


// ============================================================
// CONSTRUCTOR
// ============================================================

Camera::Camera(
    const glm::vec3& target,
    float distance,
    float fov,
    float nearPlane,
    float farPlane
)
    : position(0.0f, 0.0f, distance),
      target(target),
      distance(distance),
      yaw(0.0f),
      pitch(0.0f),
      fov(fov),
      nearPlane(nearPlane),
      farPlane(farPlane),
      aspectRatio(16.0f / 9.0f)
{
    updatePosition();
}


// ============================================================
// MOVE FORWARD / BACKWARD
// ============================================================

void Camera::moveForward(float amount)
{
    glm::vec3 direction =
        glm::normalize(
            target - position
        );


    position += direction * amount;
    target += direction * amount;
}


// ============================================================
// MOVE RIGHT / LEFT
// ============================================================

void Camera::moveRight(float amount)
{
    glm::vec3 direction =
        glm::normalize(
            target - position
        );


    glm::vec3 right =
        glm::normalize(
            glm::cross(
                direction,
                glm::vec3(0.0f, 1.0f, 0.0f)
            )
        );


    position += right * amount;
    target += right * amount;
}


// ============================================================
// ORBIT
// ============================================================

void Camera::orbit(
    float deltaX,
    float deltaY
)
{
    constexpr float sensitivity = 0.005f;


    yaw -= deltaX * sensitivity;
    pitch -= deltaY * sensitivity;


    constexpr float limit =
        glm::half_pi<float>() - 0.01f;


    pitch =
        std::clamp(
            pitch,
            -limit,
            limit
        );


    updatePosition();
}


// ============================================================
// ZOOM
// ============================================================

void Camera::zoom(float delta)
{
    constexpr float zoomSpeed = 0.5f;


    distance -=
        delta * zoomSpeed;


    distance =
        std::clamp(
            distance,
            0.5f,
            1000.0f
        );


    updatePosition();
}


// ============================================================
// UPDATE ORBIT POSITION
// ============================================================

void Camera::updatePosition()
{
    const float cosPitch =
        std::cos(pitch);

    const float sinPitch =
        std::sin(pitch);

    const float cosYaw =
        std::cos(yaw);

    const float sinYaw =
        std::sin(yaw);


    position.x =
        target.x +
        distance *
        cosPitch *
        sinYaw;


    position.y =
        target.y +
        distance *
        sinPitch;


    position.z =
        target.z +
        distance *
        cosPitch *
        cosYaw;
}


// ============================================================
// TARGET
// ============================================================

void Camera::setTarget(
    const glm::vec3& target
)
{
    this->target = target;

    updatePosition();
}


// ============================================================
// DISTANCE
// ============================================================

void Camera::setDistance(
    float distance
)
{
    this->distance =
        std::clamp(
            distance,
            0.5f,
            1000.0f
        );


    updatePosition();
}


// ============================================================
// PROJECTION
// ============================================================

void Camera::setFOV(float fov)
{
    this->fov = fov;
}


void Camera::setNearPlane(
    float nearPlane
)
{
    this->nearPlane = nearPlane;
}


void Camera::setFarPlane(
    float farPlane
)
{
    this->farPlane = farPlane;
}


void Camera::setAspectRatio(
    float aspectRatio
)
{
    this->aspectRatio = aspectRatio;
}


// ============================================================
// VIEW MATRIX
// ============================================================

glm::mat4 Camera::getViewMatrix() const
{
    return glm::lookAt(
        position,
        target,
        glm::vec3(0.0f, 1.0f, 0.0f)
    );
}


// ============================================================
// PROJECTION MATRIX
// ============================================================

glm::mat4 Camera::getProjectionMatrix() const
{
    return glm::perspective(
        glm::radians(fov),
        aspectRatio,
        nearPlane,
        farPlane
    );
}


// ============================================================
// GETTERS
// ============================================================

const glm::vec3& Camera::getPosition() const
{
    return position;
}


const glm::vec3& Camera::getTarget() const
{
    return target;
}


float Camera::getDistance() const
{
    return distance;
}


float Camera::getFOV() const
{
    return fov;
}


float Camera::getNearPlane() const
{
    return nearPlane;
}


float Camera::getFarPlane() const
{
    return farPlane;
}


float Camera::getAspectRatio() const
{
    return aspectRatio;
}