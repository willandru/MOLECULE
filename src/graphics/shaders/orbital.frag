#version 330 core

in vec3 fragmentNormal;
in vec3 fragmentPosition;

uniform vec3 orbitalColor;

out vec4 fragmentColor;

void main()
{
    vec3 normal =
        normalize(
            fragmentNormal
        );


    vec3 lightDirection =
        normalize(
            vec3(
                0.4,
                0.7,
                1.0
            )
        );


    float diffuse =
        max(
            dot(
                normal,
                lightDirection
            ),
            0.0
        );


    float ambient =
        0.30;


    float illumination =
        ambient +
        0.70 * diffuse;


    vec3 color =
        orbitalColor *
        illumination;


    fragmentColor =
        vec4(
            color,
            0.88
        );
}