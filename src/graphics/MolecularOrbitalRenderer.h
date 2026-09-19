#pragma once

#include "MolecularOrbitalIsosurface.h"

#include <glm/glm.hpp>

#include <cstddef>
#include <vector>

class Shader;

class MolecularOrbitalRenderer
{
public:
    MolecularOrbitalRenderer();

    ~MolecularOrbitalRenderer();

    MolecularOrbitalRenderer(
        const MolecularOrbitalRenderer&
    ) = delete;

    MolecularOrbitalRenderer& operator=(
        const MolecularOrbitalRenderer&
    ) = delete;

    void initialize();

    void setSurface(
        const MolecularOrbitalIsosurface::Surface& surface
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

    std::size_t getVertexCount() const;

    std::size_t getIndexCount() const;

private:
    struct GPUVertex
    {
        float x;
        float y;
        float z;

        float nx;
        float ny;
        float nz;

        float phase;
    };

    Shader* shader;

    unsigned int vao;
    unsigned int vbo;
    unsigned int ebo;

    std::size_t indexCount;
    std::size_t vertexCount;

    glm::mat4 modelMatrix;
    glm::mat4 viewMatrix;
    glm::mat4 projectionMatrix;

    glm::vec3 positiveColor;
    glm::vec3 negativeColor;

    bool initialized;

    void createBuffers();

    void destroyBuffers();

    void uploadSurface(
        const MolecularOrbitalIsosurface::Surface& surface
    );

    std::vector<GPUVertex> convertVertices(
        const MolecularOrbitalIsosurface::Surface& surface
    );
};