#include "NavigationViewController.h"

#include "AtomicDFT.h"
#include "AtomicOrbital3D.h"
#include "AtomicOrbitalAngular.h"
#include "AtomicOrbitalIsosurface.h"
#include "DFTConstants.h"
#include "ElectronicConfiguration.h"
#include "PBE96.h"
#include "RadialGrid.h"

#include "OrbitalRenderer.h"

#include "Grid.h"
#include "GridRenderer.h"

#include <glad/glad.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <cstddef>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>


namespace
{

// ============================================================
// WINDOW
// ============================================================

constexpr int WINDOW_WIDTH = 1280;
constexpr int WINDOW_HEIGHT = 720;


// ============================================================
// PERIODIC TABLE LAYOUT
// ============================================================
//
// The first 30 elements are positioned according to their
// actual periodic-table group and period.
//
// Groups:  1 ... 18
// Periods: 1 ... 4
//
// Only elements Z = 1 ... 30 are included.
//
// ============================================================

constexpr float PERIODIC_TABLE_SPACING_X = 10.0f;
constexpr float PERIODIC_TABLE_SPACING_Y = 12.0f;


// ============================================================
// ORBITAL GRID
// ============================================================

constexpr std::size_t ORBITAL_GRID_POINTS = 61;

constexpr double ORBITAL_EXTENT = 8.0;

constexpr double ORBITAL_ISOVALUE = 0.02;


// ============================================================
// VISUAL SCALE
// ============================================================

constexpr float ORBITAL_SCALE = 0.75f;


// ============================================================
// ATOM PLACEMENT
// ============================================================

struct AtomPlacement
{
    int atomicNumber;

    int period;

    int group;
};


// ============================================================
// FIRST 30 ELEMENTS
// ============================================================
//
// Period 1:
// H                                               He
//
// Period 2:
// Li  Be                              B   C   N   O   F   Ne
//
// Period 3:
// Na  Mg                              Al  Si  P   S   Cl  Ar
//
// Period 4:
// K   Ca  Sc  Ti  V   Cr  Mn  Fe  Co  Ni  Cu  Zn
//
// ============================================================

const std::vector<AtomPlacement> ATOM_PLACEMENTS =
{
    // ========================================================
    // PERIOD 1
    // ========================================================

    {  1, 1,  1 },     // H
    {  2, 1, 18 },     // He


    // ========================================================
    // PERIOD 2
    // ========================================================

    {  3, 2,  1 },     // Li
    {  4, 2,  2 },     // Be

    {  5, 2, 13 },     // B
    {  6, 2, 14 },     // C
    {  7, 2, 15 },     // N
    {  8, 2, 16 },     // O
    {  9, 2, 17 },     // F
    { 10, 2, 18 },     // Ne


    // ========================================================
    // PERIOD 3
    // ========================================================

    { 11, 3,  1 },     // Na
    { 12, 3,  2 },     // Mg

    { 13, 3, 13 },     // Al
    { 14, 3, 14 },     // Si
    { 15, 3, 15 },     // P
    { 16, 3, 16 },     // S
    { 17, 3, 17 },     // Cl
    { 18, 3, 18 },     // Ar


    // ========================================================
    // PERIOD 4
    // ========================================================

    { 19, 4,  1 },     // K
    { 20, 4,  2 },     // Ca
    { 21, 4,  3 },     // Sc
    { 22, 4,  4 },     // Ti
    { 23, 4,  5 },     // V
    { 24, 4,  6 },     // Cr
    { 25, 4,  7 },     // Mn
    { 26, 4,  8 },     // Fe
    { 27, 4,  9 },     // Co
    { 28, 4, 10 },     // Ni
    { 29, 4, 11 },     // Cu
    { 30, 4, 12 }      // Zn
};


// ============================================================
// CALCULATE ATOM POSITION
// ============================================================
//
// The complete 18-column table is centered around x = 0.
//
// Group 1  -> left
// Group 18 -> right
//
// Period 1  -> top
// Period 4  -> bottom
//
// ============================================================

glm::vec3 calculateAtomPosition(
    const AtomPlacement& placement
)
{
    constexpr float CENTER_GROUP = 9.5f;
    constexpr float CENTER_PERIOD = 2.5f;


    const float x =
        (
            static_cast<float>(placement.group) -
            CENTER_GROUP
        )
        *
        PERIODIC_TABLE_SPACING_X;


    const float y =
        (
            CENTER_PERIOD -
            static_cast<float>(placement.period)
        )
        *
        PERIODIC_TABLE_SPACING_Y;


    return glm::vec3(
        x,
        y,
        0.0f
    );
}


// ============================================================
// VISUAL ORBITAL
// ============================================================

struct VisualOrbital
{
    AtomicOrbitalIsosurface::Result surface;

    glm::vec3 position =
        glm::vec3(0.0f);

    std::string atomSymbol;

    std::string orbitalName;
};


// ============================================================
// BUILD ATOMIC ORBITALS
// ============================================================

std::vector<VisualOrbital> buildAtomOrbitals(
    int Z,
    const RadialGrid& radialGrid,
    const XCFunctional& functional,
    const glm::vec3& atomPosition
)
{
    const AtomicConfiguration configuration =
        getAtomicConfiguration(Z);


    const AtomicResult atom =
        solveAtom(
            radialGrid,
            configuration,
            functional
        );


    if (!atom.scf.converged)
    {
        throw std::runtime_error(
            "SCF no convergio para " +
            atom.symbol
        );
    }


    std::cout
        << "Atom: "
        << atom.symbol
        << " | Z = "
        << atom.Z
        << " | Electrons = "
        << atom.electrons
        << " | Converged = "
        << (
            atom.scf.converged
                ? "Yes"
                : "No"
        )
        << '\n';


    AtomicOrbital3D orbital3D;

    AtomicOrbitalIsosurface isosurface;

    std::vector<VisualOrbital> result;


    for (const AtomicOrbital& orbital :
         atom.scf.orbitals)
    {
        const int angularCount =
            AtomicOrbitalAngular::orbitalCount(
                orbital.l
            );


        for (int angularIndex = 0;
             angularIndex < angularCount;
             ++angularIndex)
        {
            const AtomicOrbitalAngularType angularType =
                AtomicOrbitalAngular::typeFromQuantumNumbers(
                    orbital.l,
                    angularIndex
                );


            const std::string angularName =
                AtomicOrbitalAngular::name(
                    angularType
                );


            const AtomicOrbital3D::Grid grid =
                orbital3D.sample(
                    radialGrid,
                    orbital,
                    angularType,
                    ORBITAL_GRID_POINTS,
                    ORBITAL_EXTENT
                );


            AtomicOrbitalIsosurface::Result surface =
                isosurface.generate(
                    grid,
                    ORBITAL_ISOVALUE
                );


            if (
                surface.positive.vertices.empty() &&
                surface.negative.vertices.empty()
            )
            {
                continue;
            }


            VisualOrbital visual;

            visual.surface =
                std::move(surface);

            visual.position =
                atomPosition;

            visual.atomSymbol =
                atom.symbol;

            visual.orbitalName =
                std::to_string(orbital.n) +
                angularName;


            std::cout
                << "  "
                << visual.atomSymbol
                << " "
                << visual.orbitalName
                << " | + vertices = "
                << visual.surface
                       .positive
                       .vertices
                       .size()
                << " | + triangles = "
                << visual.surface
                       .positive
                       .indices
                       .size() / 3
                << " | - vertices = "
                << visual.surface
                       .negative
                       .vertices
                       .size()
                << " | - triangles = "
                << visual.surface
                       .negative
                       .indices
                       .size() / 3
                << '\n';


            result.push_back(
                std::move(visual)
            );
        }
    }


    return result;
}

}


// ============================================================
// MAIN
// ============================================================

int main()
{
    try
    {
        // ====================================================
        // NAVIGATION
        // ====================================================

        NavigationViewController navigation(
            WINDOW_WIDTH,
            WINDOW_HEIGHT,
            "Molecule"
        );


        // ====================================================
        // OPENGL
        // ====================================================

        glEnable(GL_DEPTH_TEST);


        // ====================================================
        // GRID
        // ====================================================

        Grid grid;

        GridRenderer gridRenderer;

        gridRenderer.initialize(
            grid
        );


        // ====================================================
        // DFT
        // ====================================================

        PBE96 pbe96;


        RadialGrid radialGrid(
            DFTConstants::GRID_POINTS,
            DFTConstants::RMAX
        );


        // ====================================================
        // VISUAL ORBITALS
        // ====================================================

        std::vector<VisualOrbital>
            visualOrbitals;


        std::vector<
            std::unique_ptr<OrbitalRenderer>
        > orbitalRenderers;


        // ====================================================
        // CALCULATE FIRST 30 ELEMENTS
        // ====================================================

        for (const AtomPlacement& placement :
             ATOM_PLACEMENTS)
        {
            const glm::vec3 atomPosition =
                calculateAtomPosition(
                    placement
                );


            std::vector<VisualOrbital>
                atomOrbitals =
                buildAtomOrbitals(
                    placement.atomicNumber,
                    radialGrid,
                    pbe96,
                    atomPosition
                );


            for (VisualOrbital& orbital :
                 atomOrbitals)
            {
                visualOrbitals.push_back(
                    std::move(orbital)
                );
            }
        }


        // ====================================================
        // CREATE GPU RENDERERS
        // ====================================================

        orbitalRenderers.reserve(
            visualOrbitals.size()
        );


        for (const VisualOrbital& orbital :
             visualOrbitals)
        {
            std::unique_ptr<OrbitalRenderer>
                renderer =
                std::make_unique<OrbitalRenderer>();


            renderer->initialize();


            renderer->setSurface(
                orbital.surface
            );


            renderer->setPositiveColor(
                glm::vec3(
                    0.10f,
                    0.35f,
                    1.00f
                )
            );


            renderer->setNegativeColor(
                glm::vec3(
                    1.00f,
                    0.15f,
                    0.15f
                )
            );


            orbitalRenderers.push_back(
                std::move(renderer)
            );
        }


        std::cout
            << "\nTotal visual orbitals: "
            << visualOrbitals.size()
            << '\n';


        // ====================================================
        // MAIN LOOP
        // ====================================================

        while (!navigation.shouldClose())
        {
            navigation.update();


            // =================================================
            // CLEAR
            // =================================================

            glClear(
                GL_COLOR_BUFFER_BIT |
                GL_DEPTH_BUFFER_BIT
            );


            // =================================================
            // MATRICES
            // =================================================

            const glm::mat4& view =
                navigation.getViewMatrix();


            const glm::mat4& projection =
                navigation.getProjectionMatrix();


            // =================================================
            // GRID
            // =================================================

            gridRenderer.render(
                view,
                projection
            );


            // =================================================
            // ORBITALS
            // =================================================

            for (std::size_t i = 0;
                 i < orbitalRenderers.size();
                 ++i)
            {
                const VisualOrbital& orbital =
                    visualOrbitals[i];


                const glm::mat4 model =
                    glm::translate(
                        glm::mat4(1.0f),
                        orbital.position
                    )
                    *
                    glm::scale(
                        glm::mat4(1.0f),
                        glm::vec3(
                            ORBITAL_SCALE
                        )
                    );


                orbitalRenderers[i]->setModelMatrix(
                    model
                );


                orbitalRenderers[i]->setViewMatrix(
                    view
                );


                orbitalRenderers[i]->setProjectionMatrix(
                    projection
                );


                orbitalRenderers[i]->render();
            }


            // =================================================
            // PRESENT
            // =================================================

            navigation.present();
        }


        return 0;
    }
    catch (const std::exception& exception)
    {
        std::cerr
            << "\nERROR:\n"
            << exception.what()
            << '\n';

        return -1;
    }
}