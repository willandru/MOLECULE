#pragma once

#include "AtomicOrbitalIsosurface.h"
#include "OrbitalRenderer.h"

#include <glm/glm.hpp>

#include <memory>
#include <string>
#include <vector>

class RadialGrid;
class XCFunctional;

class AtomBuilder
{
public:

    AtomBuilder(
        const RadialGrid& radialGrid,
        const XCFunctional& functional,
        int atomCount,
        float spacingX,
        float spacingY
    );

    void render(
        const glm::mat4& view,
        const glm::mat4& projection,
        float scale
    );

    std::size_t orbitalCount() const;

private:

    struct VisualOrbital
    {
        AtomicOrbitalIsosurface::Result surface;

        glm::vec3 position =
            glm::vec3(0.0f);

        std::string atomSymbol;
        std::string orbitalName;

        std::unique_ptr<OrbitalRenderer> renderer;
    };

    struct AtomPlacement
    {
        int atomicNumber;
        int period;
        int group;
    };

    const RadialGrid& radialGrid_;
    const XCFunctional& functional_;

    int atomCount_;
    float spacingX_;
    float spacingY_;

    std::vector<VisualOrbital> visualOrbitals_;

    static constexpr std::size_t ORBITAL_GRID_POINTS = 61;
    static constexpr double ORBITAL_EXTENT = 8.0;
    static constexpr double ORBITAL_ISOVALUE = 0.02;

    static const std::vector<AtomPlacement>
        ATOM_PLACEMENTS;

    glm::vec3 calculateAtomPosition(
        const AtomPlacement& placement
    ) const;

    std::vector<VisualOrbital> buildAtomOrbitals(
        int atomicNumber,
        const glm::vec3& atomPosition
    ) const;

    static void initializeRenderer(
        VisualOrbital& orbital
    );
};