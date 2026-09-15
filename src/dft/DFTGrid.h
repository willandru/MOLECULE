#pragma once

#include <cstddef>
#include <glm/glm.hpp>

class DFTGrid
{
public:

    DFTGrid();

    DFTGrid(
        int nx,
        int ny,
        int nz,
        double spacing,
        const glm::dvec3& origin
    );

    void initialize(
        int nx,
        int ny,
        int nz,
        double spacing,
        const glm::dvec3& origin
    );

    // Dimensiones
    int getNx() const;
    int getNy() const;
    int getNz() const;

    // Separación entre puntos de la malla
    double getSpacing() const;

    // Posición del origen
    const glm::dvec3& getOrigin() const;

    // Número total de puntos
    std::size_t getPointCount() const;

    // Volumen elemental dV = dx * dy * dz
    double getVolumeElement() const;

    // Volumen asociado a la integración numérica
    double getVolume() const;

    // Conversión (x,y,z) -> índice lineal
    std::size_t getIndex(
        int x,
        int y,
        int z
    ) const;

    // Conversión índice lineal -> coordenada de malla
    void getCoordinates(
        std::size_t index,
        int& x,
        int& y,
        int& z
    ) const;

    // Posición física del punto
    glm::dvec3 getPosition(
        int x,
        int y,
        int z
    ) const;

    glm::dvec3 getPosition(
        std::size_t index
    ) const;

    // Comprobación de coordenadas
    bool isInside(
        int x,
        int y,
        int z
    ) const;

private:

    int nx;
    int ny;
    int nz;

    double spacing;

    glm::dvec3 origin;

    std::size_t pointCount;
};