#pragma once

#include "AtomicOrbitalIsosurface.h"

#include <glm/mat4x4.hpp>
#include <glm/vec3.hpp>

#include <cstddef>
#include <vector>

class Shader;

class OrbitalRenderer
{
public:

    OrbitalRenderer();

    ~OrbitalRenderer();

    OrbitalRenderer(
        const OrbitalRenderer&
    ) = delete;

    OrbitalRenderer& operator=(
        const OrbitalRenderer&
    ) = delete;

    OrbitalRenderer(
        OrbitalRenderer&&
    ) = delete;

    OrbitalRenderer& operator=(
        OrbitalRenderer&&
    ) = delete;


    void initialize();


    void setSurface(
        const AtomicOrbitalIsosurface::Result& surface
    );


    void setModelMatrix(
        const glm::mat4& model
    );


    void setViewMatrix(
        const glm::mat4& view
    );


    void setProjectionMatrix(
        const glm::mat4& projection
    );


    void setPositiveColor(
        const glm::vec3& color
    );


    void setNegativeColor(
        const glm::vec3& color
    );


    void render();


    bool isInitialized() const;


    std::size_t getPositiveVertexCount() const;

    std::size_t getPositiveIndexCount() const;

    std::size_t getNegativeVertexCount() const;

    std::size_t getNegativeIndexCount() const;


private:

    struct GPUVertex
    {
        float x = 0.0f;
        float y = 0.0f;
        float z = 0.0f;

        float nx = 0.0f;
        float ny = 0.0f;
        float nz = 0.0f;
    };


private:

    void createBuffers();

    void destroyBuffers();

    void uploadSurface(
        const AtomicOrbitalIsosurface::Surface& surface,
        unsigned int& vao,
        unsigned int& vbo,
        unsigned int& ebo
    );


    void renderSurface(
        const AtomicOrbitalIsosurface::Surface& surface,
        unsigned int vao,
        const glm::vec3& color
    );


    static std::vector<GPUVertex> convertVertices(
        const AtomicOrbitalIsosurface::Surface& surface
    );


private:

    Shader* shader = nullptr;

    bool initialized = false;


    unsigned int positiveVAO = 0;

    unsigned int positiveVBO = 0;

    unsigned int positiveEBO = 0;


    unsigned int negativeVAO = 0;

    unsigned int negativeVBO = 0;

    unsigned int negativeEBO = 0;


    std::size_t positiveIndexCount = 0;

    std::size_t negativeIndexCount = 0;


    glm::mat4 modelMatrix;

    glm::mat4 viewMatrix;

    glm::mat4 projectionMatrix;


    glm::vec3 positiveColor =
        glm::vec3(
            0.15f,
            0.35f,
            1.00f
        );


    glm::vec3 negativeColor =
        glm::vec3(
            1.00f,
            0.20f,
            0.20f
        );
};