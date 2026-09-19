#pragma once

#include "CartesianGrid.h"

#include <glm/glm.hpp>

#include <cstddef>
#include <vector>

class MolecularOrbitalIsosurface
{
public:
    struct Vertex
    {
        glm::vec3 position;
        glm::vec3 normal;
        float phase = 0.0f;
    };

    struct Surface
    {
        std::vector<Vertex> vertices;
        std::vector<unsigned int> indices;

        bool empty() const
        {
            return vertices.empty() || indices.empty();
        }

        void clear()
        {
            vertices.clear();
            indices.clear();
        }
    };

    static Surface generate(
        const CartesianGrid& grid,
        const std::vector<double>& psi,
        double isovalue
    );

private:
    struct GridPoint
    {
        glm::vec3 position;
        double value = 0.0;
    };

    static GridPoint makeGridPoint(
        const CartesianGrid& grid,
        const std::vector<double>& psi,
        std::size_t i,
        std::size_t j,
        std::size_t k
    );

    static glm::vec3 interpolatePosition(
        const GridPoint& a,
        const GridPoint& b,
        double isovalue
    );

    static void addTriangle(
        Surface& surface,
        const glm::vec3& a,
        const glm::vec3& b,
        const glm::vec3& c,
        float phase
    );

    static void polygonizeTetrahedron(
        Surface& surface,
        const GridPoint* points,
        double isovalue
    );
};