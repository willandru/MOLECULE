#version 330 core

in vec3 fragmentPosition;
in vec3 fragmentNormal;
in float fragmentPhase;

uniform vec3 positiveColor;
uniform vec3 negativeColor;

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

    vec3 baseColor;

    if (fragmentPhase >= 0.0)
    {
        baseColor =
            positiveColor;
    }
    else
    {
        baseColor =
            negativeColor;
    }

    float lighting =
        0.30 +
        0.70 * diffuse;

    vec3 finalColor =
        baseColor * lighting;

    FragColor =
        vec4(
            finalColor,
            1.0
        );
}