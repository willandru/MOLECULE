#include "MolecularNucleusRenderer.h"

#include "Shader.h"

#include <glad/glad.h>

#include <glm/gtc/constants.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <cmath>
#include <cstddef>
#include <stdexcept>


namespace
{

constexpr GLint POSITION_ATTRIBUTE = 0;
constexpr GLint NORMAL_ATTRIBUTE = 1;

constexpr int SPHERE_SEGMENTS = 20;
constexpr int SPHERE_RINGS = 12;

constexpr double HYDROGEN_RADIUS = 0.18;
constexpr double DEFAULT_RADIUS = 0.24;

}


MolecularNucleusRenderer::MolecularNucleusRenderer()
    : shader(nullptr),
      vao(0),
      vbo(0),
      vertexCount(0),
      modelMatrix(1.0f),
      viewMatrix(1.0f),
      projectionMatrix(1.0f),
      color(
          0.85f,
          0.85f,
          0.85f
      ),
      initialized(false)
{
}


MolecularNucleusRenderer::~MolecularNucleusRenderer()
{
    destroyBuffers();

    delete shader;
    shader = nullptr;
}


void MolecularNucleusRenderer::initialize()
{
    if (initialized)
    {
        return;
    }

    shader =
        new Shader(
            "../src/graphics/shaders/nucleus.vert",
            "../src/graphics/shaders/nucleus.frag"
        );

    createBuffers();

    initialized = true;
}


void MolecularNucleusRenderer::setMolecule(
    const Molecule& molecule
)
{
    if (!initialized)
    {
        throw std::runtime_error(
            "MolecularNucleusRenderer no ha sido inicializado."
        );
    }

    vertices.clear();

    const std::vector<Molecule::Nucleus>& nuclei =
        molecule.getNuclei();

    for (
        const Molecule::Nucleus& nucleus :
        nuclei
    )
    {
        const double radius =
            nucleus.atomicNumber == 1
                ? HYDROGEN_RADIUS
                : DEFAULT_RADIUS;

        generateSphere(
            nucleus,
            radius,
            SPHERE_SEGMENTS,
            SPHERE_RINGS
        );
    }

    uploadVertices();
}


void MolecularNucleusRenderer::setModelMatrix(
    const glm::mat4& model
)
{
    modelMatrix = model;
}


void MolecularNucleusRenderer::setViewMatrix(
    const glm::mat4& view
)
{
    viewMatrix = view;
}


void MolecularNucleusRenderer::setProjectionMatrix(
    const glm::mat4& projection
)
{
    projectionMatrix = projection;
}


void MolecularNucleusRenderer::setColor(
    const glm::vec3& newColor
)
{
    color = newColor;
}


void MolecularNucleusRenderer::render()
{
    if (!initialized)
    {
        return;
    }

    if (vertexCount == 0)
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
            "color"
        );

    glUniform3fv(
        colorLocation,
        1,
        glm::value_ptr(color)
    );

    glEnable(GL_DEPTH_TEST);

    glBindVertexArray(
        vao
    );

    glDrawArrays(
        GL_TRIANGLES,
        0,
        static_cast<GLsizei>(
            vertexCount
        )
    );

    glBindVertexArray(0);
}


bool MolecularNucleusRenderer::isInitialized() const
{
    return initialized;
}


std::size_t MolecularNucleusRenderer::getVertexCount() const
{
    return vertexCount;
}


void MolecularNucleusRenderer::createBuffers()
{
    glGenVertexArrays(
        1,
        &vao
    );

    glGenBuffers(
        1,
        &vbo
    );
}


void MolecularNucleusRenderer::destroyBuffers()
{
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

    vertexCount = 0;
}


void MolecularNucleusRenderer::generateSphere(
    const Molecule::Nucleus& nucleus,
    double radius,
    int segments,
    int rings
)
{
    const double pi =
        glm::pi<double>();

    const glm::vec3 center(
        static_cast<float>(
            nucleus.position[0]
        ),
        static_cast<float>(
            nucleus.position[1]
        ),
        static_cast<float>(
            nucleus.position[2]
        )
    );

    for (int ring = 0; ring < rings; ++ring)
    {
        const double theta0 =
            pi *
            static_cast<double>(ring) /
            static_cast<double>(rings);

        const double theta1 =
            pi *
            static_cast<double>(ring + 1) /
            static_cast<double>(rings);

        for (int segment = 0; segment < segments; ++segment)
        {
            const double phi0 =
                2.0 *
                pi *
                static_cast<double>(segment) /
                static_cast<double>(segments);

            const double phi1 =
                2.0 *
                pi *
                static_cast<double>(segment + 1) /
                static_cast<double>(segments);

            const glm::vec3 n00 =
                sphereNormal(
                    theta0,
                    phi0
                );

            const glm::vec3 n01 =
                sphereNormal(
                    theta0,
                    phi1
                );

            const glm::vec3 n10 =
                sphereNormal(
                    theta1,
                    phi0
                );

            const glm::vec3 n11 =
                sphereNormal(
                    theta1,
                    phi1
                );

            vertices.push_back(
                {
                    center +
                        static_cast<float>(radius) *
                        n00,
                    n00
                }
            );

            vertices.push_back(
                {
                    center +
                        static_cast<float>(radius) *
                        n10,
                    n10
                });

            vertices.push_back(
                {
                    center +
                        static_cast<float>(radius) *
                        n11,
                    n11
                }
            );

            vertices.push_back(
                {
                    center +
                        static_cast<float>(radius) *
                        n00,
                    n00
                }
            );

            vertices.push_back(
                {
                    center +
                        static_cast<float>(radius) *
                        n11,
                    n11
                });

            vertices.push_back(
                {
                    center +
                        static_cast<float>(radius) *
                        n01,
                    n01
                }
            );
        }
    }
}


glm::vec3 MolecularNucleusRenderer::sphereNormal(
    double theta,
    double phi
) const
{
    const float sinTheta =
        static_cast<float>(
            std::sin(theta)
        );

    const float cosTheta =
        static_cast<float>(
            std::cos(theta)
        );

    const float cosPhi =
        static_cast<float>(
            std::cos(phi)
        );

    const float sinPhi =
        static_cast<float>(
            std::sin(phi)
        );

    return glm::vec3(
        sinTheta * cosPhi,
        cosTheta,
        sinTheta * sinPhi
    );
}


void MolecularNucleusRenderer::uploadVertices()
{
    vertexCount =
        vertices.size();

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
                position
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
                normal
            )
        )
    );

    glBindVertexArray(0);
}