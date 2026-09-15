#include "HydrogenRenderer.h"
#include "GeometryConstants.h"

#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>
#include <memory>

#include <glm/geometric.hpp>
#include <glm/gtc/matrix_transform.hpp>


namespace
{
    // ========================================================
    // MATHEMATICAL CONSTANTS
    // ========================================================

    constexpr double PI =
        3.14159265358979323846;


    // ========================================================
    // ISOSURFACE GEOMETRY
    // ========================================================

    constexpr int LATITUDE_SEGMENTS =
        96;

    constexpr int LONGITUDE_SEGMENTS =
        192;


    // ========================================================
    // SHADERS
    // ========================================================

    const char* VERTEX_SHADER =
        "../src/graphics/shaders/hydrogen.vert";


    const char* FRAGMENT_SHADER =
        "../src/graphics/shaders/hydrogen.frag";
}


// ============================================================
// CONSTRUCTOR
// ============================================================

HydrogenRenderer::HydrogenRenderer()
    : vao(0),
      vbo(0),
      ebo(0),
      shader(nullptr),
      position(0.0f),
      scale(1.0f),
      isovalue(0.02),
      maximumDensity(0.0),
      geometryRadius(0.0),
      initialized(false),
      geometryBuilt(false)
{
}


// ============================================================
// DESTRUCTOR
// ============================================================

HydrogenRenderer::~HydrogenRenderer()
{
    destroyOpenGLObjects();
}


// ============================================================
// INITIALIZE
// ============================================================

bool HydrogenRenderer::initialize()
{
    if (initialized)
    {
        return true;
    }


    try
    {
        shader =
            std::make_unique<Shader>(
                VERTEX_SHADER,
                FRAGMENT_SHADER
            );
    }
    catch (const std::exception& exception)
    {
        std::cerr
            << "ERROR: No se pudo crear el shader "
            << "de HydrogenRenderer.\n"
            << exception.what()
            << std::endl;

        shader.reset();

        return false;
    }


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


    if (vao == 0 ||
        vbo == 0 ||
        ebo == 0)
    {
        std::cerr
            << "ERROR: No se pudieron crear "
            << "los buffers de HydrogenRenderer."
            << std::endl;

        destroyOpenGLObjects();

        return false;
    }


    initialized = true;

    return true;
}


// ============================================================
// BUILD GEOMETRY
// ============================================================

void HydrogenRenderer::buildGeometry(
    const HydrogenAtom& hydrogen
)
{
    if (!initialized)
    {
        std::cerr
            << "ERROR: HydrogenRenderer no esta "
            << "inicializado."
            << std::endl;

        return;
    }


    if (!hydrogen.isDFTConverged())
    {
        std::cerr
            << "ERROR: El DFT del hidrogeno "
            << "no convergio."
            << std::endl;

        geometryBuilt = false;

        return;
    }


    generateIsosurface(
        hydrogen
    );


    if (vertices.empty() ||
        indices.empty())
    {
        std::cerr
            << "ERROR: No se genero geometria "
            << "para el hidrogeno."
            << std::endl;

        geometryBuilt = false;

        return;
    }


    // ========================================================
    // VERTEX ARRAY
    // ========================================================

    glBindVertexArray(
        vao
    );


    // ========================================================
    // VERTEX BUFFER
    // ========================================================

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


    // ========================================================
    // INDEX BUFFER
    // ========================================================

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


    // ========================================================
    // POSITION ATTRIBUTE
    // ========================================================

    glEnableVertexAttribArray(0);


    glVertexAttribPointer(
        0,
        3,
        GL_FLOAT,
        GL_FALSE,
        sizeof(Vertex),
        reinterpret_cast<void*>(
            offsetof(
                Vertex,
                position
            )
        )
    );


    // ========================================================
    // NORMAL ATTRIBUTE
    // ========================================================

    glEnableVertexAttribArray(1);


    glVertexAttribPointer(
        1,
        3,
        GL_FLOAT,
        GL_FALSE,
        sizeof(Vertex),
        reinterpret_cast<void*>(
            offsetof(
                Vertex,
                normal
            )
        )
    );


    glBindVertexArray(0);


    geometryBuilt = true;
}


// ============================================================
// RENDER
// ============================================================

void HydrogenRenderer::render(
    const glm::mat4& view,
    const glm::mat4& projection
)
{
    if (!initialized ||
        !geometryBuilt ||
        !shader)
    {
        return;
    }


    glm::mat4 model =
        glm::translate(
            glm::mat4(1.0f),
            position
        );


    model =
        glm::scale(
            model,
            glm::vec3(scale)
        );


    shader->use();


    shader->setMat4(
        "uModel",
        model
    );


    shader->setMat4(
        "uView",
        view
    );


    shader->setMat4(
        "uProjection",
        projection
    );


    glBindVertexArray(
        vao
    );


    glDrawElements(
        GL_TRIANGLES,
        static_cast<GLsizei>(
            indices.size()
        ),
        GL_UNSIGNED_INT,
        nullptr
    );


    glBindVertexArray(0);
}


// ============================================================
// POSITION
// ============================================================

void HydrogenRenderer::setPosition(
    const glm::vec3& position
)
{
    this->position =
        position;
}


// ============================================================
// SCALE
// ============================================================

void HydrogenRenderer::setScale(
    float scale
)
{
    if (scale <= 0.0f)
    {
        scale = 1.0f;
    }


    this->scale =
        scale;
}


// ============================================================
// ISOVALUE
// ============================================================

void HydrogenRenderer::setIsovalue(
    double value
)
{
    if (value <= 0.0)
    {
        value = 0.02;
    }


    if (value > 1.0)
    {
        value = 1.0;
    }


    isovalue =
        value;
}


double HydrogenRenderer::getIsovalue() const
{
    return isovalue;
}


// ============================================================
// INFORMATION
// ============================================================

std::size_t HydrogenRenderer::getVertexCount() const
{
    return vertices.size();
}


std::size_t HydrogenRenderer::getTriangleCount() const
{
    return indices.size() / 3;
}


// ============================================================
// GENERATE ISOSURFACE
// ============================================================

void HydrogenRenderer::generateIsosurface(
    const HydrogenAtom& hydrogen
)
{
    vertices.clear();

    indices.clear();


    maximumDensity = 0.0;

    geometryRadius = 0.0;


    const auto& radialGrid =
        hydrogen.getRadialGrid();


    const auto& density =
        hydrogen.getDensity();


    if (radialGrid.empty() ||
        density.empty() ||
        radialGrid.size() != density.size())
    {
        return;
    }


    // ========================================================
    // MAXIMUM DENSITY
    // ========================================================

    for (double value : density)
    {
        if (value > maximumDensity)
        {
            maximumDensity =
                value;
        }
    }


    if (maximumDensity <= 0.0)
    {
        return;
    }


    // ========================================================
    // ISOSURFACE RADIUS
    // ========================================================

    //
    // geometryRadius remains in Bohr here.
    //

    geometryRadius =
        findIsosurfaceRadius(
            radialGrid,
            density
        );


    if (geometryRadius <= 0.0)
    {
        return;
    }


    generateVertices(
        geometryRadius
    );


    generateIndices();

    calculateNormals();
}


// ============================================================
// SAMPLE DENSITY
// ============================================================

double HydrogenRenderer::sampleDensity(
    double radius,
    const std::vector<double>& radialGrid,
    const std::vector<double>& density
) const
{
    if (radialGrid.empty() ||
        density.empty())
    {
        return 0.0;
    }


    if (radius <= radialGrid.front())
    {
        return density.front();
    }


    if (radius >= radialGrid.back())
    {
        return density.back();
    }


    for (std::size_t i = 1;
         i < radialGrid.size();
         ++i)
    {
        if (radius <= radialGrid[i])
        {
            const double r0 =
                radialGrid[i - 1];


            const double r1 =
                radialGrid[i];


            const double d0 =
                density[i - 1];


            const double d1 =
                density[i];


            const double t =
                (radius - r0) /
                (r1 - r0);


            return d0 +
                   t * (d1 - d0);
        }
    }


    return density.back();
}


// ============================================================
// FIND ISOSURFACE RADIUS
// ============================================================

double HydrogenRenderer::findIsosurfaceRadius(
    const std::vector<double>& radialGrid,
    const std::vector<double>& density
) const
{
    if (radialGrid.empty() ||
        density.empty())
    {
        return 0.0;
    }


    const double target =
        maximumDensity *
        isovalue;


    for (std::size_t i = 1;
         i < radialGrid.size();
         ++i)
    {
        const double d0 =
            density[i - 1];


        const double d1 =
            density[i];


        if (d0 >= target &&
            d1 <= target)
        {
            const double r0 =
                radialGrid[i - 1];


            const double r1 =
                radialGrid[i];


            const double denominator =
                d1 - d0;


            if (std::abs(denominator) <
                1.0e-15)
            {
                return r1;
            }


            const double t =
                (target - d0) /
                denominator;


            return r0 +
                   t * (r1 - r0);
        }
    }


    return radialGrid.back();
}


// ============================================================
// GENERATE VERTICES
// ============================================================

void HydrogenRenderer::generateVertices(
    double radius
)
{
    vertices.clear();


    // ========================================================
    // BOHR -> ANGSTROM
    // ========================================================

    //
    // DFT:
    //
    //     radius = Bohr
    //
    // Geometry:
    //
    //     1 internal unit = 1 Å
    //
    // Therefore:
    //
    //     radius(Å) =
    //     radius(Bohr) *
    //     0.529177210903
    //

    const double visualRadius =
        radius *
        GeometryConstants::BOHR_TO_ANGSTROM;


    const int columns =
        LONGITUDE_SEGMENTS + 1;


    vertices.reserve(
        static_cast<std::size_t>(
            (LATITUDE_SEGMENTS + 1) *
            columns
        )
    );


    for (int latitude = 0;
         latitude <= LATITUDE_SEGMENTS;
         ++latitude)
    {
        const double v =
            static_cast<double>(latitude) /
            static_cast<double>(
                LATITUDE_SEGMENTS
            );


        const double theta =
            v * PI;


        const double sinTheta =
            std::sin(theta);


        const double cosTheta =
            std::cos(theta);


        for (int longitude = 0;
             longitude <= LONGITUDE_SEGMENTS;
             ++longitude)
        {
            const double u =
                static_cast<double>(longitude) /
                static_cast<double>(
                    LONGITUDE_SEGMENTS
                );


            const double phi =
                u * 2.0 * PI;


            const double sinPhi =
                std::sin(phi);


            const double cosPhi =
                std::cos(phi);


            glm::vec3 normal(
                static_cast<float>(
                    sinTheta * cosPhi
                ),

                static_cast<float>(
                    cosTheta
                ),

                static_cast<float>(
                    sinTheta * sinPhi
                )
            );


            normal =
                glm::normalize(
                    normal
                );


            Vertex vertex;


            vertex.position =
                normal *
                static_cast<float>(
                    visualRadius
                );


            vertex.normal =
                normal;


            vertices.push_back(
                vertex
            );
        }
    }
}


// ============================================================
// GENERATE INDICES
// ============================================================

void HydrogenRenderer::generateIndices()
{
    indices.clear();


    const int columns =
        LONGITUDE_SEGMENTS + 1;


    indices.reserve(
        static_cast<std::size_t>(
            LATITUDE_SEGMENTS *
            LONGITUDE_SEGMENTS *
            6
        )
    );


    for (int latitude = 0;
         latitude < LATITUDE_SEGMENTS;
         ++latitude)
    {
        for (int longitude = 0;
             longitude < LONGITUDE_SEGMENTS;
             ++longitude)
        {
            const unsigned int current =
                static_cast<unsigned int>(
                    latitude * columns +
                    longitude
                );


            const unsigned int next =
                current + 1;


            const unsigned int below =
                static_cast<unsigned int>(
                    (latitude + 1) *
                    columns +
                    longitude
                );


            const unsigned int belowNext =
                below + 1;


            indices.push_back(
                current
            );

            indices.push_back(
                below
            );

            indices.push_back(
                next
            );


            indices.push_back(
                next
            );

            indices.push_back(
                below
            );

            indices.push_back(
                belowNext
            );
        }
    }
}


// ============================================================
// CALCULATE NORMALS
// ============================================================

void HydrogenRenderer::calculateNormals()
{
    for (Vertex& vertex : vertices)
    {
        if (glm::length(vertex.position) > 0.0f)
        {
            vertex.normal =
                glm::normalize(
                    vertex.position
                );
        }
        else
        {
            vertex.normal =
                glm::vec3(
                    0.0f,
                    1.0f,
                    0.0f
                );
        }
    }
}


// ============================================================
// CLEANUP
// ============================================================

void HydrogenRenderer::destroyOpenGLObjects()
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


    shader.reset();


    initialized = false;

    geometryBuilt = false;
}