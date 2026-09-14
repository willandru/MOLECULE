#include "Shader.h"

#include <glm/gtc/type_ptr.hpp>

#include <fstream>
#include <sstream>
#include <stdexcept>


// ============================================================
// CONSTRUCTOR
// ============================================================

Shader::Shader(
    const std::string& vertexPath,
    const std::string& fragmentPath
)
    : program(0)
{
    const std::string vertexSource =
        readFile(vertexPath);

    const std::string fragmentSource =
        readFile(fragmentPath);


    GLuint vertexShader =
        compileShader(
            GL_VERTEX_SHADER,
            vertexSource
        );


    GLuint fragmentShader =
        compileShader(
            GL_FRAGMENT_SHADER,
            fragmentSource
        );


    program =
        createProgram(
            vertexShader,
            fragmentShader
        );


    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);
}


// ============================================================
// DESTRUCTOR
// ============================================================

Shader::~Shader()
{
    if (program != 0)
    {
        glDeleteProgram(program);

        program = 0;
    }
}


// ============================================================
// USE
// ============================================================

void Shader::use() const
{
    glUseProgram(program);
}


// ============================================================
// SET MATRIX
// ============================================================

void Shader::setMat4(
    const std::string& name,
    const glm::mat4& matrix
) const
{
    GLint location =
        glGetUniformLocation(
            program,
            name.c_str()
        );


    if (location >= 0)
    {
        glUniformMatrix4fv(
            location,
            1,
            GL_FALSE,
            glm::value_ptr(matrix)
        );
    }
}


// ============================================================
// GET PROGRAM
// ============================================================

GLuint Shader::getProgram() const
{
    return program;
}


// ============================================================
// READ FILE
// ============================================================

std::string Shader::readFile(
    const std::string& path
)
{
    std::ifstream file(path);


    if (!file.is_open())
    {
        throw std::runtime_error(
            "Failed to open shader file: " +
            path
        );
    }


    std::stringstream buffer;

    buffer << file.rdbuf();


    return buffer.str();
}


// ============================================================
// COMPILE SHADER
// ============================================================

GLuint Shader::compileShader(
    GLenum type,
    const std::string& source
)
{
    GLuint shader =
        glCreateShader(type);


    const char* sourceCode =
        source.c_str();


    glShaderSource(
        shader,
        1,
        &sourceCode,
        nullptr
    );


    glCompileShader(shader);


    GLint success = GL_FALSE;


    glGetShaderiv(
        shader,
        GL_COMPILE_STATUS,
        &success
    );


    if (!success)
    {
        GLint logLength = 0;


        glGetShaderiv(
            shader,
            GL_INFO_LOG_LENGTH,
            &logLength
        );


        std::string log(
            logLength,
            '\0'
        );


        glGetShaderInfoLog(
            shader,
            logLength,
            nullptr,
            log.data()
        );


        glDeleteShader(shader);


        throw std::runtime_error(
            "Shader compilation failed:\n" +
            log
        );
    }


    return shader;
}


// ============================================================
// CREATE PROGRAM
// ============================================================

GLuint Shader::createProgram(
    GLuint vertexShader,
    GLuint fragmentShader
)
{
    GLuint program =
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


    GLint success = GL_FALSE;


    glGetProgramiv(
        program,
        GL_LINK_STATUS,
        &success
    );


    if (!success)
    {
        GLint logLength = 0;


        glGetProgramiv(
            program,
            GL_INFO_LOG_LENGTH,
            &logLength
        );


        std::string log(
            logLength,
            '\0'
        );


        glGetProgramInfoLog(
            program,
            logLength,
            nullptr,
            log.data()
        );


        glDeleteProgram(program);


        throw std::runtime_error(
            "Shader linking failed:\n" +
            log
        );
    }


    return program;
}
