#include "MolecularOrbitalRenderer.h"

#include "Shader.h"

#include <glad/glad.h>

#include <glm/gtc/type_ptr.hpp>

#include <stdexcept>
#include <vector>


namespace
{

constexpr GLint POSITION_ATTRIBUTE = 0;
constexpr GLint NORMAL_ATTRIBUTE = 1;
constexpr GLint PHASE_ATTRIBUTE = 2;

}


MolecularOrbitalRenderer::MolecularOrbitalRenderer()
    : shader(nullptr),
      vao(0),
      vbo(0),
      ebo(0),
      indexCount(0),
      vertexCount(0),
      modelMatrix(1.0f),
      viewMatrix(1.0f),
      projectionMatrix(1.0f),
      positiveColor(
          0.15f,
          0.35f,
          1.0f
      ),
      negativeColor(
          1.0f,
          0.20f,
          0.20f
      ),
      initialized(false)
{
}


MolecularOrbitalRenderer::~MolecularOrbitalRenderer()
{
    destroyBuffers();

    delete shader;
    shader = nullptr;
}


void MolecularOrbitalRenderer::initialize()
{
    if (initialized)
    {
        return;
    }

    shader = new Shader(
        "../src/graphics/shaders/molecular.vert",
        "../src/graphics/shaders/molecular.frag"
    );

    createBuffers();

    initialized = true;
}


void MolecularOrbitalRenderer::setSurface(
    const MolecularOrbitalIsosurface::Surface& surface
)
{
    if (!initialized)
    {
        throw std::runtime_error(
            "MolecularOrbitalRenderer no ha sido inicializado."
        );
    }

    uploadSurface(surface);

    vertexCount =
        surface.vertices.size();

    indexCount =
        surface.indices.size();
}


void MolecularOrbitalRenderer::setModelMatrix(
    const glm::mat4& model
)
{
    modelMatrix = model;
}


void MolecularOrbitalRenderer::setViewMatrix(
    const glm::mat4& view
)
{
    viewMatrix = view;
}


void MolecularOrbitalRenderer::setProjectionMatrix(
    const glm::mat4& projection
)
{
    projectionMatrix = projection;
}


void MolecularOrbitalRenderer::setPositiveColor(
    const glm::vec3& color
)
{
    positiveColor = color;
}


void MolecularOrbitalRenderer::setNegativeColor(
    const glm::vec3& color
)
{
    negativeColor = color;
}


void MolecularOrbitalRenderer::render()
{
    if (!initialized)
    {
        return;
    }

    if (indexCount == 0)
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

    const GLint positiveColorLocation =
        glGetUniformLocation(
            program,
            "positiveColor"
        );

    const GLint negativeColorLocation =
        glGetUniformLocation(
            program,
            "negativeColor"
        );

    glUniform3fv(
        positiveColorLocation,
        1,
        glm::value_ptr(
            positiveColor
        )
    );

    glUniform3fv(
        negativeColorLocation,
        1,
        glm::value_ptr(
            negativeColor
        )
    );

    glEnable(GL_DEPTH_TEST);

    glDisable(GL_CULL_FACE);

    glBindVertexArray(
        vao
    );

    glDrawElements(
        GL_TRIANGLES,
        static_cast<GLsizei>(
            indexCount
        ),
        GL_UNSIGNED_INT,
        nullptr
    );

    glBindVertexArray(0);

    glEnable(GL_CULL_FACE);
}


bool MolecularOrbitalRenderer::isInitialized() const
{
    return initialized;
}


std::size_t MolecularOrbitalRenderer::getVertexCount() const
{
    return vertexCount;
}


std::size_t MolecularOrbitalRenderer::getIndexCount() const
{
    return indexCount;
}


void MolecularOrbitalRenderer::createBuffers()
{
    glGenVertexArrays(
        1,
        &vao
    );

    glGenBuffers(
        1,
        &vbo
    );

    glGenBuffers(
        1,
        &ebo
    );
}


void MolecularOrbitalRenderer::destroyBuffers()
{
    if (ebo != 0)
    {
        glDeleteBuffers(
            1,
            &ebo
        );

        ebo = 0;
    }

    if (vbo != 0)
    {
        glDeleteBuffers(
            1,
            &vbo
        );

        vbo = 0;
    }

    if (vao != 0)
    {
        glDeleteVertexArrays(
            1,
            &vao
        );

        vao = 0;
    }

    indexCount = 0;
    vertexCount = 0;
}


void MolecularOrbitalRenderer::uploadSurface(
    const MolecularOrbitalIsosurface::Surface& surface
)
{
    const std::vector<GPUVertex> vertices =
        convertVertices(surface);

    glBindVertexArray(
        vao
    );

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
            offsetof(
                GPUVertex,
                x
            )
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
            offsetof(
                GPUVertex,
                nx
            )
        )
    );

    glEnableVertexAttribArray(
        PHASE_ATTRIBUTE
    );

    glVertexAttribPointer(
        PHASE_ATTRIBUTE,
        1,
        GL_FLOAT,
        GL_FALSE,
        sizeof(GPUVertex),
        reinterpret_cast<void*>(
            offsetof(
                GPUVertex,
                phase
            )
        )
    );

    glBindVertexArray(0);
}


std::vector<MolecularOrbitalRenderer::GPUVertex>
MolecularOrbitalRenderer::convertVertices(
    const MolecularOrbitalIsosurface::Surface& surface
)
{
    std::vector<GPUVertex> vertices;

    vertices.reserve(
        surface.vertices.size()
    );

    for (
        const MolecularOrbitalIsosurface::Vertex& vertex :
        surface.vertices
    )
    {
        GPUVertex gpuVertex;

        gpuVertex.x =
            vertex.position.x;

        gpuVertex.y =
            vertex.position.y;

        gpuVertex.z =
            vertex.position.z;

        gpuVertex.nx =
            vertex.normal.x;

        gpuVertex.ny =
            vertex.normal.y;

        gpuVertex.nz =
            vertex.normal.z;

        gpuVertex.phase =
            vertex.phase;

        vertices.push_back(
            gpuVertex
        );
    }

    return vertices;
}