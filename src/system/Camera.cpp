#include "Camera.h"

#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>

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
      pitch(glm::radians(18.0f)),
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
//
// Movimiento sobre el plano horizontal de la cámara.
//
// W  -> adelante
// S  -> atrás
//
// El movimiento NO modifica distance.
// Por lo tanto W/S nunca actúa como zoom.
//

void Camera::moveForward(
    float amount
)
{
    glm::vec3 forward =
        target - position;


    // Movimiento horizontal.
    forward.y = 0.0f;


    const float length =
        glm::length(forward);


    if (length < 1.0e-6f)
    {
        return;
    }


    forward /= length;


    position +=
        forward * amount;

    target +=
        forward * amount;
}


// ============================================================
// MOVE RIGHT / LEFT
// ============================================================
//
// A -> izquierda
// D -> derecha
//

void Camera::moveRight(
    float amount
)
{
    glm::vec3 forward =
        target - position;


    forward.y = 0.0f;


    const float length =
        glm::length(forward);


    if (length < 1.0e-6f)
    {
        return;
    }


    forward /= length;


    glm::vec3 right =
        glm::cross(
            forward,
            glm::vec3(
                0.0f,
                1.0f,
                0.0f
            )
        );


    const float rightLength =
        glm::length(right);


    if (rightLength < 1.0e-6f)
    {
        return;
    }


    right /= rightLength;


    position +=
        right * amount;

    target +=
        right * amount;
}


// ============================================================
// MOVE VERTICAL
// ============================================================
//
// Q -> abajo
// E -> arriba
//

void Camera::moveVertical(
    float amount
)
{
    position.y += amount;
    target.y += amount;
}


// ============================================================
// ORBIT
// ============================================================

void Camera::orbit(
    float deltaX,
    float deltaY
)
{
    constexpr float sensitivity =
        0.005f;


    yaw -=
        deltaX * sensitivity;

    pitch -=
        deltaY * sensitivity;


    constexpr float limit =
        glm::half_pi<float>() - 0.05f;


    pitch =
        std::clamp(
            pitch,
            -limit,
            limit
        );


    updatePosition();
}


// ============================================================
// SET ORIENTATION
// ============================================================

void Camera::setOrientation(
    float yaw,
    float pitch
)
{
    this->yaw =
        yaw;


    constexpr float limit =
        glm::half_pi<float>() - 0.05f;


    this->pitch =
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
//
// La rueda modifica exclusivamente la distancia orbital.
//
// Esto es independiente de W/S.
//

void Camera::zoom(
    float delta
)
{
    constexpr float zoomSpeed =
        1.0f;


    distance -=
        delta * zoomSpeed;


    distance =
        std::clamp(
            distance,
            1.0f,
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
        target.y -
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
    this->target =
        target;


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
            1.0f,
            1000.0f
        );


    updatePosition();
}


// ============================================================
// FOV
// ============================================================

void Camera::setFOV(
    float fov
)
{
    this->fov =
        fov;
}


// ============================================================
// NEAR PLANE
// ============================================================

void Camera::setNearPlane(
    float nearPlane
)
{
    this->nearPlane =
        nearPlane;
}


// ============================================================
// FAR PLANE
// ============================================================

void Camera::setFarPlane(
    float farPlane
)
{
    this->farPlane =
        farPlane;
}


// ============================================================
// ASPECT RATIO
// ============================================================

void Camera::setAspectRatio(
    float aspectRatio
)
{
    this->aspectRatio =
        aspectRatio;
}


// ============================================================
// VIEW MATRIX
// ============================================================

glm::mat4 Camera::getViewMatrix() const
{
    return glm::lookAt(
        position,
        target,
        glm::vec3(
            0.0f,
            1.0f,
            0.0f
        )
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
// POSITION
// ============================================================

const glm::vec3& Camera::getPosition() const
{
    return position;
}


// ============================================================
// TARGET
// ============================================================

const glm::vec3& Camera::getTarget() const
{
    return target;
}


// ============================================================
// DISTANCE
// ============================================================

float Camera::getDistance() const
{
    return distance;
}


// ============================================================
// FOV
// ============================================================

float Camera::getFOV() const
{
    return fov;
}


// ============================================================
// NEAR PLANE
// ============================================================

float Camera::getNearPlane() const
{
    return nearPlane;
}


// ============================================================
// FAR PLANE
// ============================================================

float Camera::getFarPlane() const
{
    return farPlane;
}


// ============================================================
// ASPECT RATIO
// ============================================================

float Camera::getAspectRatio() const
{
    return aspectRatio;
}