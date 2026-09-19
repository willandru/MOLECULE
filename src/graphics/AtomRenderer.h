#pragma once

#include <glm/glm.hpp>

#include "AtomicDFT.h"
#include "RadialGrid.h"

class Shader;

class AtomRenderer
{
public:
    AtomRenderer();
    ~AtomRenderer();

    AtomRenderer(const AtomRenderer&) = delete;
    AtomRenderer& operator=(const AtomRenderer&) = delete;

    void initialize(
        const RadialGrid& grid,
        const AtomicResult& atom,
        const glm::vec3& position,
        float densityThreshold
    );

    void render(
        const glm::mat4& view,
        const glm::mat4& projection
    );

private:
    void buildGeometry(
        const RadialGrid& grid,
        const AtomicResult& atom,
        const glm::vec3& position,
        float densityThreshold
    );

private:
    unsigned int vao;
    unsigned int vbo;
    unsigned int ebo;

    int indexCount;

    Shader* shader;
};