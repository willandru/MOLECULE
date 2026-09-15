#version 330 core

in vec3 vNormal;
in vec3 vWorldPosition;

out vec4 FragColor;

void main()
{
    vec3 normal =
        normalize(vNormal);

    vec3 lightDirection =
        normalize(
            vec3(
                0.6,
                0.8,
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
        0.25;

    float intensity =
        ambient +
        0.75 * diffuse;

    vec3 baseColor =
        vec3(
            0.85,
            0.90,
            1.00
        );

    vec3 color =
        baseColor *
        intensity;

    FragColor =
        vec4(
            color,
            1.0
        );
}