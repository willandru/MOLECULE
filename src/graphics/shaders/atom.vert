#version 330 core

layout(location = 0) in vec3 aPosition;
layout(location = 1) in vec3 aNormal;


// ============================================================
// MATRICES
// ============================================================

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;


// ============================================================
// OUTPUT
// ============================================================

out vec3 WorldPosition;
out vec3 WorldNormal;


// ============================================================
// MAIN
// ============================================================

void main()
{
    vec4 worldPosition =
        model * vec4(aPosition, 1.0);

    WorldPosition =
        worldPosition.xyz;

    WorldNormal =
        normalize(
            mat3(transpose(inverse(model))) *
            aNormal
        );

    gl_Position =
        projection *
        view *
        worldPosition;
}