#include "OrbitalRenderer.h"

#include "Shader.h"

#include <glad/glad.h>

#include <glm/gtc/type_ptr.hpp>

#include <stdexcept>
#include <vector>


namespace
{

constexpr GLint POSITION_ATTRIBUTE = 0;
constexpr GLint NORMAL_ATTRIBUTE = 1;

}


OrbitalRenderer::OrbitalRenderer()
    : modelMatrix(1.0f),
      viewMatrix(1.0f),
      projectionMatrix(1.0f)
{
}


OrbitalRenderer::~OrbitalRenderer()
{
    destroyBuffers();

    delete shader;
    shader = nullptr;
}


void OrbitalRenderer::initialize()
{
    if (initialized)
    {
        return;
    }

    shader = new Shader(
        "../src/graphics/shaders/orbital.vert",
        "../src/graphics/shaders/orbital.frag"
    );

    createBuffers();

    initialized = true;
}


void OrbitalRenderer::setSurface(
    const AtomicOrbitalIsosurface::Result& surface
)
{
    if (!initialized)
    {
        throw std::runtime_error(
            "OrbitalRenderer no ha sido inicializado."
        );
    }

    uploadSurface(
        surface.positive,
        positiveVAO,
        positiveVBO,
        positiveEBO
    );

    positiveIndexCount =
        surface.positive.indices.size();

    uploadSurface(
        surface.negative,
        negativeVAO,
        negativeVBO,
        negativeEBO
    );

    negativeIndexCount =
        surface.negative.indices.size();
}


void OrbitalRenderer::setModelMatrix(
    const glm::mat4& model
)
{
    modelMatrix = model;
}


void OrbitalRenderer::setViewMatrix(
    const glm::mat4& view
)
{
    viewMatrix = view;
}


void OrbitalRenderer::setProjectionMatrix(
    const glm::mat4& projection
)
{
    projectionMatrix = projection;
}


void OrbitalRenderer::setPositiveColor(
    const glm::vec3& color
)
{
    positiveColor = color;
}


void OrbitalRenderer::setNegativeColor(
    const glm::vec3& color
)
{
    negativeColor = color;
}


void OrbitalRenderer::render()
{
    if (!initialized)
    {
        return;
    }

    if (positiveIndexCount == 0 &&
        negativeIndexCount == 0)
    {
        return;
    }

    shader->use();

    shader->setMat4(
        "model",
        modelMatrix
    );

    shader->setMat4(
        "view",
        viewMatrix
    );

    shader->setMat4(
        "projection",
        projectionMatrix
    );

    const GLuint program =
        shader->getProgram();

    const GLint colorLocation =
        glGetUniformLocation(
            program,
            "orbitalColor"
        );

    glEnable(GL_DEPTH_TEST);

    glDisable(GL_CULL_FACE);

    if (positiveIndexCount > 0)
    {
        glUniform3fv(
            colorLocation,
            1,
            glm::value_ptr(
                positiveColor
            )
        );

        glBindVertexArray(
            positiveVAO
        );

        glDrawElements(
            GL_TRIANGLES,
            static_cast<GLsizei>(
                positiveIndexCount
            ),
            GL_UNSIGNED_INT,
            nullptr
        );
    }

    if (negativeIndexCount > 0)
    {
        glUniform3fv(
            colorLocation,
            1,
            glm::value_ptr(
                negativeColor
            )
        );

        glBindVertexArray(
            negativeVAO
        );

        glDrawElements(
            GL_TRIANGLES,
            static_cast<GLsizei>(
                negativeIndexCount
            ),
            GL_UNSIGNED_INT,
            nullptr
        );
    }

    glBindVertexArray(0);

    glEnable(GL_CULL_FACE);
}


bool OrbitalRenderer::isInitialized() const
{
    return initialized;
}


std::size_t OrbitalRenderer::getPositiveVertexCount() const
{
    return 0;
}


std::size_t OrbitalRenderer::getPositiveIndexCount() const
{
    return positiveIndexCount;
}


std::size_t OrbitalRenderer::getNegativeVertexCount() const
{
    return 0;
}


std::size_t OrbitalRenderer::getNegativeIndexCount() const
{
    return negativeIndexCount;
}


void OrbitalRenderer::createBuffers()
{
    glGenVertexArrays(
        1,
        &positiveVAO
    );

    glGenBuffers(
        1,
        &positiveVBO
    );

    glGenBuffers(
        1,
        &positiveEBO
    );

    glGenVertexArrays(
        1,
        &negativeVAO
    );

    glGenBuffers(
        1,
        &negativeVBO
    );

    glGenBuffers(
        1,
        &negativeEBO
    );
}


void OrbitalRenderer::destroyBuffers()
{
    if (positiveEBO != 0)
    {
        glDeleteBuffers(
            1,
            &positiveEBO
        );

        positiveEBO = 0;
    }

    if (positiveVBO != 0)
    {
        glDeleteBuffers(
            1,
            &positiveVBO
        );

        positiveVBO = 0;
    }

    if (positiveVAO != 0)
    {
        glDeleteVertexArrays(
            1,
            &positiveVAO
        );

        positiveVAO = 0;
    }

    if (negativeEBO != 0)
    {
        glDeleteBuffers(
            1,
            &negativeEBO
        );

        negativeEBO = 0;
    }

    if (negativeVBO != 0)
    {
        glDeleteBuffers(
            1,
            &negativeVBO
        );

        negativeVBO = 0;
    }

    if (negativeVAO != 0)
    {
        glDeleteVertexArrays(
            1,
            &negativeVAO
        );

        negativeVAO = 0;
    }

    positiveIndexCount = 0;
    negativeIndexCount = 0;
}


void OrbitalRenderer::uploadSurface(
    const AtomicOrbitalIsosurface::Surface& surface,
    unsigned int& vao,
    unsigned int& vbo,
    unsigned int& ebo
)
{
    const std::vector<GPUVertex> vertices =
        convertVertices(surface);

    glBindVertexArray(vao);

    glBindBuffer(
        GL_ARRAY_BUFFER,
        vbo
    );

    glBufferData(
        GL_ARRAY_BUFFER,
        static_cast<GLsizeiptr>(
            vertices.size() *
            sizeof(GPUVertex)
        ),
        vertices.empty()
            ? nullptr
            : vertices.data(),
        GL_STATIC_DRAW
    );

    glBindBuffer(
        GL_ELEMENT_ARRAY_BUFFER,
        ebo
    );

    glBufferData(
        GL_ELEMENT_ARRAY_BUFFER,
        static_cast<GLsizeiptr>(
            surface.indices.size() *
            sizeof(unsigned int)
        ),
        surface.indices.empty()
            ? nullptr
            : surface.indices.data(),
        GL_STATIC_DRAW
    );

    glEnableVertexAttribArray(
        POSITION_ATTRIBUTE
    );

    glVertexAttribPointer(
        POSITION_ATTRIBUTE,
        3,
        GL_FLOAT,
        GL_FALSE,
        sizeof(GPUVertex),
        reinterpret_cast<void*>(
            0
        )
    );

    glEnableVertexAttribArray(
        NORMAL_ATTRIBUTE
    );

    glVertexAttribPointer(
        NORMAL_ATTRIBUTE,
        3,
        GL_FLOAT,
        GL_FALSE,
        sizeof(GPUVertex),
        reinterpret_cast<void*>(
            3 * sizeof(float)
        )
    );

    glBindVertexArray(0);
}


void OrbitalRenderer::renderSurface(
    const AtomicOrbitalIsosurface::Surface&,
    unsigned int,
    const glm::vec3&
)
{
}


std::vector<OrbitalRenderer::GPUVertex>
OrbitalRenderer::convertVertices(
    const AtomicOrbitalIsosurface::Surface& surface
)
{
    std::vector<GPUVertex> vertices;

    vertices.reserve(
        surface.vertices.size()
    );

    for (
        const AtomicOrbitalIsosurface::Vertex& vertex :
        surface.vertices
    )
    {
        GPUVertex gpuVertex;

        gpuVertex.x =
            static_cast<float>(
                vertex.x
            );

        gpuVertex.y =
            static_cast<float>(
                vertex.y
            );

        gpuVertex.z =
            static_cast<float>(
                vertex.z
            );

        gpuVertex.nx =
            static_cast<float>(
                vertex.nx
            );

        gpuVertex.ny =
            static_cast<float>(
                vertex.ny
            );

        gpuVertex.nz =
            static_cast<float>(
                vertex.nz
            );

        vertices.push_back(
            gpuVertex
        );
    }

    return vertices;
}