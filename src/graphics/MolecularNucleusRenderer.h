#pragma once

#include "Molecule.h"

#include <glm/glm.hpp>

#include <cstddef>
#include <vector>

class Shader;

class MolecularNucleusRenderer
{
public:
    MolecularNucleusRenderer();

    ~MolecularNucleusRenderer();

    MolecularNucleusRenderer(
        const MolecularNucleusRenderer&
    ) = delete;

    MolecularNucleusRenderer& operator=(
        const MolecularNucleusRenderer&
    ) = delete;

    void initialize();

    void setMolecule(
        const Molecule& molecule
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

    void setColor(
        const glm::vec3& color
    );

    void render();

    bool isInitialized() const;

    std::size_t getVertexCount() const;

private:
    struct GPUVertex
    {
        glm::vec3 position;
        glm::vec3 normal;
    };

    Shader* shader;

    unsigned int vao;
    unsigned int vbo;

    std::size_t vertexCount;

    glm::mat4 modelMatrix;
    glm::mat4 viewMatrix;
    glm::mat4 projectionMatrix;

    glm::vec3 color;

    bool initialized;

    std::vector<GPUVertex> vertices;

    void createBuffers();

    void destroyBuffers();

    void generateSphere(
        const Molecule::Nucleus& nucleus,
        double radius,
        int segments,
        int rings
    );

    glm::vec3 sphereNormal(
        double theta,
        double phi
    ) const;

    void uploadVertices();
};