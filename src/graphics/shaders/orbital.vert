#version 330 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

out vec3 fragmentNormal;
out vec3 fragmentPosition;

void main()
{
    vec4 worldPosition =
        model *
        vec4(position, 1.0);

    fragmentPosition =
        worldPosition.xyz;

    fragmentNormal =
        mat3(
            transpose(
                inverse(model)
            )
        )
        *
        normal;

    gl_Position =
        projection *
        view *
        worldPosition;
}