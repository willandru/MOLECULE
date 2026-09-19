#include "AtomRenderer.h"

#include "Shader.h"

#include <glad/glad.h>
#include <glm/gtc/constants.hpp>

#include <cmath>
#include <stdexcept>
#include <vector>

namespace
{

struct Vertex
{
    glm::vec3 position;
};

float findDensityRadius(
    const RadialGrid& grid,
    const AtomicResult& atom,
    float densityThreshold
)
{
    const std::vector<double>& r =
        grid.coordinates();

    const std::vector<double>& density =
        atom.scf.density;

    if (r.size() != density.size())
    {
        throw std::invalid_argument(
            "La malla radial y la densidad del atomo deben tener el mismo tamano."
        );
    }

    if (r.size() < 2)
    {
        throw std::invalid_argument(
            "La malla radial debe contener al menos dos puntos."
        );
    }

    if (densityThreshold <= 0.0f)
    {
        throw std::invalid_argument(
            "El umbral de densidad debe ser mayor que cero."
        );
    }

    for (std::size_t i = 0; i < r.size(); ++i)
    {
        if (density[i] <= static_cast<double>(densityThreshold))
        {
            if (i == 0)
            {
                return static_cast<float>(r[i]);
            }

            const double r0 = r[i - 1];
            const double r1 = r[i];

            const double rho0 = density[i - 1];
            const double rho1 = density[i];

            if (std::abs(rho1 - rho0) < 1.0e-30)
            {
                return static_cast<float>(r1);
            }

            const double t =
                (
                    static_cast<double>(densityThreshold) -
                    rho0
                ) /
                (
                    rho1 - rho0
                );

            const double radius =
                r0 +
                t * (r1 - r0);

            return static_cast<float>(radius);
        }
    }

    return static_cast<float>(r.back());
}

}

AtomRenderer::AtomRenderer()
    : vao(0),
      vbo(0),
      ebo(0),
      indexCount(0),
      shader(nullptr)
{
}

AtomRenderer::~AtomRenderer()
{
    if (ebo != 0)
    {
        glDeleteBuffers(
            1,
            &ebo
        );
    }

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

    delete shader;
}

void AtomRenderer::initialize(
    const RadialGrid& grid,
    const AtomicResult& atom,
    const glm::vec3& position,
    float densityThreshold
)
{
    if (!atom.scf.converged)
    {
        throw std::runtime_error(
            "No se puede renderizar un atomo cuyo calculo DFT no convergio."
        );
    }

    buildGeometry(
        grid,
        atom,
        position,
        densityThreshold
    );

    shader =
        new Shader(
            "../src/graphics/shaders/atom.vert",
            "../src/graphics/shaders/atom.frag"
        );
}

void AtomRenderer::buildGeometry(
    const RadialGrid& grid,
    const AtomicResult& atom,
    const glm::vec3& position,
    float densityThreshold
)
{
    const float radius =
        findDensityRadius(
            grid,
            atom,
            densityThreshold
        );

    constexpr int latitudeSegments = 32;
    constexpr int longitudeSegments = 64;

    std::vector<Vertex> vertices;
    std::vector<unsigned int> indices;

    vertices.reserve(
        static_cast<std::size_t>(
            (latitudeSegments + 1) *
            (longitudeSegments + 1)
        )
    );

    indices.reserve(
        static_cast<std::size_t>(
            latitudeSegments *
            longitudeSegments *
            6
        )
    );

    const float pi =
        glm::pi<float>();

    for (int latitude = 0;
         latitude <= latitudeSegments;
         ++latitude)
    {
        const float theta =
            pi *
            static_cast<float>(latitude) /
            static_cast<float>(latitudeSegments);

        const float sinTheta =
            std::sin(theta);

        const float cosTheta =
            std::cos(theta);

        for (int longitude = 0;
             longitude <= longitudeSegments;
             ++longitude)
        {
            const float phi =
                2.0f *
                pi *
                static_cast<float>(longitude) /
                static_cast<float>(longitudeSegments);

            const float sinPhi =
                std::sin(phi);

            const float cosPhi =
                std::cos(phi);

            glm::vec3 vertexPosition;

            vertexPosition.x =
                position.x +
                radius *
                sinTheta *
                cosPhi;

            vertexPosition.y =
                position.y +
                radius *
                cosTheta;

            vertexPosition.z =
                position.z +
                radius *
                sinTheta *
                sinPhi;

            vertices.push_back({
                vertexPosition
            });
        }
    }

    const unsigned int stride =
        static_cast<unsigned int>(
            longitudeSegments + 1
        );

    for (int latitude = 0;
         latitude < latitudeSegments;
         ++latitude)
    {
        for (int longitude = 0;
             longitude < longitudeSegments;
             ++longitude)
        {
            const unsigned int current =
                static_cast<unsigned int>(
                    latitude
                ) *
                stride +
                static_cast<unsigned int>(
                    longitude
                );

            const unsigned int next =
                current +
                stride;

            indices.push_back(current);
            indices.push_back(next);
            indices.push_back(current + 1);

            indices.push_back(current + 1);
            indices.push_back(next);
            indices.push_back(next + 1);
        }
    }

    indexCount =
        static_cast<int>(
            indices.size()
        );

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
            sizeof(Vertex)
        ),
        vertices.data(),
        GL_STATIC_DRAW
    );

    glBindBuffer(
        GL_ELEMENT_ARRAY_BUFFER,
        ebo
    );

    glBufferData(
        GL_ELEMENT_ARRAY_BUFFER,
        static_cast<GLsizeiptr>(
            indices.size() *
            sizeof(unsigned int)
        ),
        indices.data(),
        GL_STATIC_DRAW
    );

    glEnableVertexAttribArray(0);

    glVertexAttribPointer(
        0,
        3,
        GL_FLOAT,
        GL_FALSE,
        sizeof(Vertex),
        nullptr
    );

    glBindVertexArray(0);
}

void AtomRenderer::render(
    const glm::mat4& view,
    const glm::mat4& projection
)
{
    if (shader == nullptr ||
        vao == 0 ||
        indexCount <= 0)
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

    glBindVertexArray(
        vao
    );

    glDrawElements(
        GL_TRIANGLES,
        indexCount,
        GL_UNSIGNED_INT,
        nullptr
    );

    glBindVertexArray(0);
}