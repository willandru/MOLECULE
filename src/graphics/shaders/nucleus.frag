#version 330 core

in vec3 fragmentNormal;

uniform vec3 color;

out vec4 FragColor;

void main()
{
    vec3 normal =
        normalize(fragmentNormal);

    vec3 lightDirection =
        normalize(
            vec3(
                0.4,
                0.7,
                0.6
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

    float lighting =
        0.30 +
        0.70 * diffuse;

    FragColor =
        vec4(
            color * lighting,
            1.0
        );
}