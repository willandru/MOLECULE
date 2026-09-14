#pragma once

#include <glad/glad.h>
#include <glm/glm.hpp>


class Grid;
class Shader;


class GridRenderer
{
public:

    GridRenderer();

    ~GridRenderer();

    GridRenderer(const GridRenderer&) = delete;
    GridRenderer& operator=(const GridRenderer&) = delete;


    void initialize(
        const Grid& grid
    );


    void render(
        const glm::mat4& view,
        const glm::mat4& projection
    );


private:

    GLuint vao;
    GLuint vbo;

    GLsizei vertexCount;

    Shader* shader;
};