#pragma once

#include <string>

#include <glad/glad.h>
#include <glm/glm.hpp>


class Shader
{
public:

    Shader(
        const std::string& vertexPath,
        const std::string& fragmentPath
    );

    ~Shader();

    Shader(const Shader&) = delete;
    Shader& operator=(const Shader&) = delete;


    void use() const;


    void setMat4(
        const std::string& name,
        const glm::mat4& matrix
    ) const;


    GLuint getProgram() const;


private:

    GLuint program;


    static std::string readFile(
        const std::string& path
    );


    static GLuint compileShader(
        GLenum type,
        const std::string& source
    );


    static GLuint createProgram(
        GLuint vertexShader,
        GLuint fragmentShader
    );
};
