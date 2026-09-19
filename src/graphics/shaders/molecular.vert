#version 330 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in float phase;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

out vec3 fragmentPosition;
out vec3 fragmentNormal;
out float fragmentPhase;

void main()
{
    vec4 worldPosition =
        model * vec4(position, 1.0);

    fragmentPosition =
        worldPosition.xyz;

    fragmentNormal =
        mat3(transpose(inverse(model))) *
        normal;

    fragmentPhase =
        phase;

    gl_Position =
        projection *
        view *
        worldPosition;
}