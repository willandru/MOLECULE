#include "GridRenderer.h"

#include "Grid.h"
#include "Shader.h"

#include <glm/glm.hpp>

#include <vector>


// ============================================================
// CONSTRUCTOR
// ============================================================

GridRenderer::GridRenderer()
    : vao(0),
      vbo(0),
      vertexCount(0),
      shader(nullptr)
{
}


// ============================================================
// DESTRUCTOR
// ============================================================

GridRenderer::~GridRenderer()
{
    delete shader;


    if (vbo != 0)
    {
        glDeleteBuffers(
            1,
            &vbo
        );
    }


    if (vao != 0)
    {
        glDeleteVertexArrays(
            1,
            &vao
        );
    }
}


// ============================================================
// INITIALIZE
// ============================================================

void GridRenderer::initialize(
    const Grid& grid
)
{
    std::vector<glm::vec3> vertices;


    const glm::vec3& min =
        grid.getMin();


    const glm::vec3& max =
        grid.getMax();


    const int divisions =
        grid.getDivisions();


    const float spacing =
        grid.getSpacing();


    // ========================================================
    // XY PLANE
    //
    // Z = 0
    // ========================================================

    for (int i = 0; i <= divisions; ++i)
    {
        const float coordinate =
            min.x +
            static_cast<float>(i) * spacing;


        // ----------------------------------------------------
        // LÍNEA PARALELA AL EJE Y
        // ----------------------------------------------------

        vertices.emplace_back(
            coordinate,
            min.y,
            0.0f
        );


        vertices.emplace_back(
            coordinate,
            max.y,
            0.0f
        );


        // ----------------------------------------------------
        // LÍNEA PARALELA AL EJE X
        // ----------------------------------------------------

        vertices.emplace_back(
            min.x,
            coordinate,
            0.0f
        );


        vertices.emplace_back(
            max.x,
            coordinate,
            0.0f
        );
    }


    vertexCount =
        static_cast<GLsizei>(
            vertices.size()
        );


    // ========================================================
    // VAO
    // ========================================================

    glGenVertexArrays(
        1,
        &vao
    );


    glBindVertexArray(vao);


    // ========================================================
    // VBO
    // ========================================================

    glGenBuffers(
        1,
        &vbo
    );


    glBindBuffer(
        GL_ARRAY_BUFFER,
        vbo
    );


    glBufferData(
        GL_ARRAY_BUFFER,
        vertices.size() * sizeof(glm::vec3),
        vertices.data(),
        GL_STATIC_DRAW
    );


    // ========================================================
    // POSITION ATTRIBUTE
    // ========================================================

    glVertexAttribPointer(
        0,
        3,
        GL_FLOAT,
        GL_FALSE,
        sizeof(glm::vec3),
        nullptr
    );


    glEnableVertexAttribArray(0);


    // ========================================================
    // UNBIND
    // ========================================================

    glBindBuffer(
        GL_ARRAY_BUFFER,
        0
    );


    glBindVertexArray(0);


    // ========================================================
    // SHADER
    // ========================================================

    shader =
        new Shader(
            "../src/Graphics/shaders/grid.vert",
            "../src/Graphics/shaders/grid.frag"
        );
}


// ============================================================
// RENDER
// ============================================================

void GridRenderer::render(
    const glm::mat4& view,
    const glm::mat4& projection
)
{
    if (shader == nullptr)
    {
        return;
    }


    shader->use();


    shader->setMat4(
        "view",
        view
    );


    shader->setMat4(
        "projection",
        projection
    );


    glBindVertexArray(vao);


    glDrawArrays(
        GL_LINES,
        0,
        vertexCount
    );


    glBindVertexArray(0);
}
