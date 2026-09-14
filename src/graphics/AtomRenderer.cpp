#include "AtomRenderer.h"
#include "AtomColorsConstants.h"

#include <glad/glad.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdio>


// ============================================================
// SPHERE VERTEX
// ============================================================

struct AtomVertex
{
    glm::vec3 position;
    glm::vec3 normal;
};


// ============================================================
// SHADER FILE
// ============================================================

static std::string loadShaderFile(
    const char* path
)
{
    std::ifstream file(path);

    if (!file.is_open())
    {
        std::printf(
            "ERROR: Could not open shader: %s\n",
            path
        );

        return {};
    }

    std::stringstream buffer;

    buffer << file.rdbuf();

    return buffer.str();
}


// ============================================================
// SHADER COMPILATION
// ============================================================

static unsigned int compileShader(
    unsigned int type,
    const char* source
)
{
    unsigned int shader =
        glCreateShader(type);


    glShaderSource(
        shader,
        1,
        &source,
        nullptr
    );


    glCompileShader(shader);


    int success = 0;

    glGetShaderiv(
        shader,
        GL_COMPILE_STATUS,
        &success
    );


    if (!success)
    {
        char infoLog[1024];

        glGetShaderInfoLog(
            shader,
            1024,
            nullptr,
            infoLog
        );

        std::printf(
            "ERROR: Atom shader compilation failed:\n%s\n",
            infoLog
        );

        glDeleteShader(shader);

        return 0;
    }


    return shader;
}


// ============================================================
// SHADER PROGRAM
// ============================================================

static unsigned int createShaderProgram(
    const char* vertexPath,
    const char* fragmentPath
)
{
    std::string vertexSource =
        loadShaderFile(vertexPath);


    std::string fragmentSource =
        loadShaderFile(fragmentPath);


    if (vertexSource.empty() ||
        fragmentSource.empty())
    {
        return 0;
    }


    unsigned int vertexShader =
        compileShader(
            GL_VERTEX_SHADER,
            vertexSource.c_str()
        );


    if (vertexShader == 0)
    {
        return 0;
    }


    unsigned int fragmentShader =
        compileShader(
            GL_FRAGMENT_SHADER,
            fragmentSource.c_str()
        );


    if (fragmentShader == 0)
    {
        glDeleteShader(vertexShader);

        return 0;
    }


    unsigned int program =
        glCreateProgram();


    glAttachShader(
        program,
        vertexShader
    );


    glAttachShader(
        program,
        fragmentShader
    );


    glLinkProgram(program);


    int success = 0;

    glGetProgramiv(
        program,
        GL_LINK_STATUS,
        &success
    );


    if (!success)
    {
        char infoLog[1024];

        glGetProgramInfoLog(
            program,
            1024,
            nullptr,
            infoLog
        );

        std::printf(
            "ERROR: Atom shader linking failed:\n%s\n",
            infoLog
        );

        glDeleteShader(vertexShader);
        glDeleteShader(fragmentShader);
        glDeleteProgram(program);

        return 0;
    }


    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);


    return program;
}


// ============================================================
// CONSTRUCTOR
// ============================================================

AtomRenderer::AtomRenderer()
    : m_VAO(0),
      m_VBO(0),
      m_EBO(0),
      m_shaderProgram(0),
      m_indexCount(0),
      m_initialized(false),
      m_lightPosition(5.0f, 8.0f, 10.0f)
{
}


// ============================================================
// DESTRUCTOR
// ============================================================

AtomRenderer::~AtomRenderer()
{
    shutdown();
}


// ============================================================
// INITIALIZE
// ============================================================

bool AtomRenderer::initialize()
{
    if (m_initialized)
    {
        return true;
    }


    // --------------------------------------------------------
    // SHADER
    // --------------------------------------------------------

    m_shaderProgram =
        createShaderProgram(
            "../src/graphics/shaders/atom.vert",
            "../src/graphics/shaders/atom.frag"
        );


    if (m_shaderProgram == 0)
    {
        return false;
    }


    // --------------------------------------------------------
    // SPHERE
    // --------------------------------------------------------

    createSphere();


    m_initialized = true;

    return true;
}


// ============================================================
// CREATE SPHERE
// ============================================================

void AtomRenderer::createSphere()
{
    constexpr int stacks = 32;
    constexpr int sectors = 32;

    constexpr float PI =
        3.14159265358979323846f;


    std::vector<AtomVertex> vertices;
    std::vector<unsigned int> indices;


    // --------------------------------------------------------
    // VERTICES
    // --------------------------------------------------------

    for (int i = 0; i <= stacks; ++i)
    {
        const float stackAngle =
            PI * 0.5f -
            static_cast<float>(i) *
            PI /
            static_cast<float>(stacks);


        const float xy =
            std::cos(stackAngle);


        const float z =
            std::sin(stackAngle);


        for (int j = 0; j <= sectors; ++j)
        {
            const float sectorAngle =
                static_cast<float>(j) *
                2.0f *
                PI /
                static_cast<float>(sectors);


            const float x =
                xy *
                std::cos(sectorAngle);


            const float y =
                xy *
                std::sin(sectorAngle);


            glm::vec3 position(
                x,
                y,
                z
            );


            glm::vec3 normal =
                glm::normalize(position);


            vertices.push_back(
                {
                    position,
                    normal
                }
            );
        }
    }


    // --------------------------------------------------------
    // INDICES
    // --------------------------------------------------------

    for (int i = 0; i < stacks; ++i)
    {
        for (int j = 0; j < sectors; ++j)
        {
            const unsigned int first =
                i *
                (sectors + 1) +
                j;


            const unsigned int second =
                first +
                sectors +
                1;


            indices.push_back(first);
            indices.push_back(second);
            indices.push_back(first + 1);


            indices.push_back(second);
            indices.push_back(second + 1);
            indices.push_back(first + 1);
        }
    }


    m_indexCount =
        static_cast<int>(
            indices.size()
        );


    // --------------------------------------------------------
    // VAO
    // --------------------------------------------------------

    glGenVertexArrays(
        1,
        &m_VAO
    );


    glGenBuffers(
        1,
        &m_VBO
    );


    glGenBuffers(
        1,
        &m_EBO
    );


    glBindVertexArray(
        m_VAO
    );


    // --------------------------------------------------------
    // VBO
    // --------------------------------------------------------

    glBindBuffer(
        GL_ARRAY_BUFFER,
        m_VBO
    );


    glBufferData(
        GL_ARRAY_BUFFER,
        static_cast<GLsizeiptr>(
            vertices.size() *
            sizeof(AtomVertex)
        ),
        vertices.data(),
        GL_STATIC_DRAW
    );


    // --------------------------------------------------------
    // EBO
    // --------------------------------------------------------

    glBindBuffer(
        GL_ELEMENT_ARRAY_BUFFER,
        m_EBO
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


    // --------------------------------------------------------
    // POSITION
    // --------------------------------------------------------

    glVertexAttribPointer(
        0,
        3,
        GL_FLOAT,
        GL_FALSE,
        sizeof(AtomVertex),
        reinterpret_cast<void*>(0)
    );


    glEnableVertexAttribArray(0);


    // --------------------------------------------------------
    // NORMAL
    // --------------------------------------------------------

    glVertexAttribPointer(
        1,
        3,
        GL_FLOAT,
        GL_FALSE,
        sizeof(AtomVertex),
        reinterpret_cast<void*>(
            sizeof(glm::vec3)
        )
    );


    glEnableVertexAttribArray(1);


    glBindVertexArray(0);
}


// ============================================================
// RENDER
// ============================================================

void AtomRenderer::render(
    const Atom& atom,
    const glm::mat4& view,
    const glm::mat4& projection
)
{
    if (!m_initialized)
    {
        return;
    }


    glUseProgram(
        m_shaderProgram
    );


    // --------------------------------------------------------
    // MODEL
    // --------------------------------------------------------

    glm::mat4 model(1.0f);


    model =
        glm::translate(
            model,
            atom.getPosition()
        );


    model =
        glm::scale(
            model,
            glm::vec3(
                atom.getCovalentRadius()
            )
        );


    // --------------------------------------------------------
    // UNIFORMS
    // --------------------------------------------------------

    const int modelLocation =
        glGetUniformLocation(
            m_shaderProgram,
            "model"
        );


    const int viewLocation =
        glGetUniformLocation(
            m_shaderProgram,
            "view"
        );


    const int projectionLocation =
        glGetUniformLocation(
            m_shaderProgram,
            "projection"
        );


    const int lightPositionLocation =
        glGetUniformLocation(
            m_shaderProgram,
            "lightPosition"
        );


    const int atomColorLocation =
        glGetUniformLocation(
            m_shaderProgram,
            "atomColor"
        );


    // --------------------------------------------------------
    // MATRICES
    // --------------------------------------------------------

    glUniformMatrix4fv(
        modelLocation,
        1,
        GL_FALSE,
        &model[0][0]
    );


    glUniformMatrix4fv(
        viewLocation,
        1,
        GL_FALSE,
        &view[0][0]
    );


    glUniformMatrix4fv(
        projectionLocation,
        1,
        GL_FALSE,
        &projection[0][0]
    );


    // --------------------------------------------------------
    // LIGHT
    // --------------------------------------------------------

    glUniform3fv(
        lightPositionLocation,
        1,
        &m_lightPosition[0]
    );


    // --------------------------------------------------------
    // COLOR
    // --------------------------------------------------------

    const glm::vec3 color =
        AtomColors::get(
            atom.getAtomicNumber()
        );


    glUniform3fv(
        atomColorLocation,
        1,
        &color[0]
    );


    // --------------------------------------------------------
    // DRAW
    // --------------------------------------------------------

    glBindVertexArray(
        m_VAO
    );


    glDrawElements(
        GL_TRIANGLES,
        m_indexCount,
        GL_UNSIGNED_INT,
        nullptr
    );


    glBindVertexArray(0);
}


// ============================================================
// SHUTDOWN
// ============================================================

void AtomRenderer::shutdown()
{
    if (m_VAO != 0)
    {
        glDeleteVertexArrays(
            1,
            &m_VAO
        );

        m_VAO = 0;
    }


    if (m_VBO != 0)
    {
        glDeleteBuffers(
            1,
            &m_VBO
        );

        m_VBO = 0;
    }


    if (m_EBO != 0)
    {
        glDeleteBuffers(
            1,
            &m_EBO
        );

        m_EBO = 0;
    }


    if (m_shaderProgram != 0)
    {
        glDeleteProgram(
            m_shaderProgram
        );

        m_shaderProgram = 0;
    }


    m_indexCount = 0;

    m_initialized = false;
}