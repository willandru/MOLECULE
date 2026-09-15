#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include <glad/glad.h>

#include <glm/mat4x4.hpp>
#include <glm/vec3.hpp>

#include "HydrogenAtom.h"
#include "Shader.h"


class HydrogenRenderer
{
public:

    HydrogenRenderer();

    ~HydrogenRenderer();


    HydrogenRenderer(
        const HydrogenRenderer&
    ) = delete;


    HydrogenRenderer& operator=(
        const HydrogenRenderer&
    ) = delete;


    // ========================================================
    // INITIALIZATION
    // ========================================================

    bool initialize();


    // ========================================================
    // GEOMETRY
    // ========================================================

    void buildGeometry(
        const HydrogenAtom& hydrogen
    );


    // ========================================================
    // RENDER
    // ========================================================

    void render(
        const glm::mat4& view,
        const glm::mat4& projection
    );


    // ========================================================
    // TRANSFORM
    // ========================================================

    void setPosition(
        const glm::vec3& position
    );

    void setScale(
        float scale
    );


    // ========================================================
    // ISOSURFACE
    // ========================================================

    void setIsovalue(
        double value
    );

    double getIsovalue() const;


    // ========================================================
    // INFORMATION
    // ========================================================

    std::size_t getVertexCount() const;

    std::size_t getTriangleCount() const;


private:

    struct Vertex
    {
        glm::vec3 position;

        glm::vec3 normal;
    };


    // ========================================================
    // GEOMETRY
    // ========================================================

    void generateIsosurface(
        const HydrogenAtom& hydrogen
    );


    double sampleDensity(
        double radius,
        const std::vector<double>& radialGrid,
        const std::vector<double>& density
    ) const;


    double findIsosurfaceRadius(
        const std::vector<double>& radialGrid,
        const std::vector<double>& density
    ) const;


    void generateVertices(
        double radius
    );


    void generateIndices();

    void calculateNormals();


    // ========================================================
    // OPENGL
    // ========================================================

    void destroyOpenGLObjects();


private:

    GLuint vao;

    GLuint vbo;

    GLuint ebo;


    std::unique_ptr<Shader> shader;


    std::vector<Vertex> vertices;

    std::vector<unsigned int> indices;


    glm::vec3 position;

    float scale;


    double isovalue;

    double maximumDensity;

    double geometryRadius;


    bool initialized;

    bool geometryBuilt;
};