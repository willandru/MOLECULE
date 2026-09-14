#pragma once

#include <glm/glm.hpp>


class Grid
{
public:
    Grid();

    float getSize() const;
    int getDivisions() const;
    float getSpacing() const;

    const glm::vec3& getMin() const;
    const glm::vec3& getMax() const;

private:
    float size;
    int divisions;
    float spacing;

    glm::vec3 min;
    glm::vec3 max;
};