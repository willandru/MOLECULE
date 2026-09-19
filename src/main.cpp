#include "NavigationViewController.h"

#include "CartesianGrid.h"
#include "DFTData.h"
#include "Grid.h"
#include "GridRenderer.h"
#include "MolecularNucleusRenderer.h"
#include "MolecularOrbitalIsosurface.h"
#include "MolecularOrbitalRenderer.h"
#include "Molecule.h"
#include "PZ81.h"
#include "SelfConsistentFieldMolecule.h"

#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <cmath>
#include <exception>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

namespace
{

constexpr int WINDOW_WIDTH = 1280;
constexpr int WINDOW_HEIGHT = 720;

constexpr double ISOVALUE = 0.04;

constexpr float H2_POSITION_X = -5.0f;
constexpr float H2_POSITION_Y = 0.0f;

constexpr float H2O_POSITION_X = 5.0f;
constexpr float H2O_POSITION_Y = 0.0f;

constexpr float POSITIVE_RED = 0.15f;
constexpr float POSITIVE_GREEN = 0.45f;
constexpr float POSITIVE_BLUE = 1.0f;

constexpr float NEGATIVE_RED = 1.0f;
constexpr float NEGATIVE_GREEN = 0.20f;
constexpr float NEGATIVE_BLUE = 0.20f;

constexpr float NUCLEUS_RED = 0.90f;
constexpr float NUCLEUS_GREEN = 0.90f;
constexpr float NUCLEUS_BLUE = 0.90f;

double maximumAbsoluteValue(
    const std::vector<double>& values
)
{
    double maximum = 0.0;

    for (const double value : values)
    {
        maximum =
            std::max(
                maximum,
                std::abs(value)
            );
    }

    return maximum;
}

}

int main()
{
    try
    {
        // =====================================================
        // H2
        // =====================================================

        Molecule h2;

        h2.addNucleus(
            1,
            -0.7,
            0.0,
            0.0
        );

        h2.addNucleus(
            1,
            0.7,
            0.0,
            0.0
        );

        CartesianGrid h2Grid(
            21,
            21,
            21,
            -3.8,
            4.2,
            -3.8,
            4.2,
            -3.8,
            4.2
        );


        // =====================================================
        // H2O
        // =====================================================

        Molecule h2o;

        h2o.addNucleus(
            8,
            0.0,
            0.0,
            0.0
        );

        h2o.addNucleus(
            1,
            1.430,
            0.0,
            1.108
        );

        h2o.addNucleus(
            1,
            -1.430,
            0.0,
            1.108
        );

        CartesianGrid h2oGrid(
            21,
            21,
            21,
            -3.8,
            4.2,
            -3.8,
            4.2,
            -3.8,
            4.2
        );


        // =====================================================
        // DFT
        // =====================================================

        PZ81 functional;

        MolecularResult h2Result =
            solveMolecularSelfConsistentField(
                h2Grid,
                h2,
                0,
                1,
                functional
            );

        MolecularResult h2oResult =
            solveMolecularSelfConsistentField(
                h2oGrid,
                h2o,
                0,
                1,
                functional
            );


        if (!h2Result.scf.converged)
        {
            throw std::runtime_error(
                "H2 molecular SCF did not converge."
            );
        }

        if (!h2oResult.scf.converged)
        {
            throw std::runtime_error(
                "H2O molecular SCF did not converge."
            );
        }


        // =====================================================
        // ORBITAL INFORMATION
        // =====================================================

        std::cout
            << "\n========================================\n"
            << "H2 MOLECULAR ORBITALS\n"
            << "========================================\n";

        for (
            std::size_t i = 0;
            i < h2Result.scf.molecularOrbitals.size();
            ++i
        )
        {
            const MolecularOrbital& orbital =
                h2Result.scf.molecularOrbitals[i];

            const double maximum =
                maximumAbsoluteValue(
                    orbital.psi
                );

            std::cout
                << "MO "
                << i
                << " | eigenvalue = "
                << orbital.eigenvalue
                << " | electrons = "
                << orbital.electrons
                << " | max|psi| = "
                << maximum
                << '\n';
        }


        std::cout
            << "\n========================================\n"
            << "H2O MOLECULAR ORBITALS\n"
            << "========================================\n";

        for (
            std::size_t i = 0;
            i < h2oResult.scf.molecularOrbitals.size();
            ++i
        )
        {
            const MolecularOrbital& orbital =
                h2oResult.scf.molecularOrbitals[i];

            const double maximum =
                maximumAbsoluteValue(
                    orbital.psi
                );

            std::cout
                << "MO "
                << i
                << " | eigenvalue = "
                << orbital.eigenvalue
                << " | electrons = "
                << orbital.electrons
                << " | max|psi| = "
                << maximum
                << '\n';
        }


        // =====================================================
        // ISOSURFACES
        // =====================================================

        std::vector<
            MolecularOrbitalIsosurface::Surface
        > h2Surfaces;

        h2Surfaces.reserve(
            h2Result.scf.molecularOrbitals.size()
        );

        for (
            const MolecularOrbital& orbital :
            h2Result.scf.molecularOrbitals
        )
        {
            MolecularOrbitalIsosurface::Surface surface =
                MolecularOrbitalIsosurface::generate(
                    h2Grid,
                    orbital.psi,
                    ISOVALUE
                );

            h2Surfaces.push_back(
                std::move(surface)
            );
        }


        std::vector<
            MolecularOrbitalIsosurface::Surface
        > h2oSurfaces;

        h2oSurfaces.reserve(
            h2oResult.scf.molecularOrbitals.size()
        );

        for (
            const MolecularOrbital& orbital :
            h2oResult.scf.molecularOrbitals
        )
        {
            MolecularOrbitalIsosurface::Surface surface =
                MolecularOrbitalIsosurface::generate(
                    h2oGrid,
                    orbital.psi,
                    ISOVALUE
                );

            h2oSurfaces.push_back(
                std::move(surface)
            );
        }


        // =====================================================
        // WINDOW
        // =====================================================

        NavigationViewController navigation(
            WINDOW_WIDTH,
            WINDOW_HEIGHT,
            "Molecular Orbitals"
        );


        // =====================================================
        // GRID
        // =====================================================

        Grid grid;

        GridRenderer gridRenderer;

        gridRenderer.initialize(
            grid
        );


        // =====================================================
        // H2 ORBITAL RENDERERS
        // =====================================================

        std::vector<
            std::unique_ptr<MolecularOrbitalRenderer>
        > h2OrbitalRenderers;

        for (
            const MolecularOrbitalIsosurface::Surface& surface :
            h2Surfaces
        )
        {
            if (surface.empty())
            {
                h2OrbitalRenderers.push_back(
                    nullptr
                );

                continue;
            }

            auto renderer =
                std::make_unique<MolecularOrbitalRenderer>();

            renderer->initialize();

            renderer->setSurface(
                surface
            );

            renderer->setPositiveColor(
                glm::vec3(
                    POSITIVE_RED,
                    POSITIVE_GREEN,
                    POSITIVE_BLUE
                )
            );

            renderer->setNegativeColor(
                glm::vec3(
                    NEGATIVE_RED,
                    NEGATIVE_GREEN,
                    NEGATIVE_BLUE
                )
            );

            renderer->setModelMatrix(
                glm::translate(
                    glm::mat4(1.0f),
                    glm::vec3(
                        H2_POSITION_X,
                        H2_POSITION_Y,
                        0.0f
                    )
                )
            );

            h2OrbitalRenderers.push_back(
                std::move(renderer)
            );
        }


        // =====================================================
        // H2O ORBITAL RENDERERS
        // =====================================================

        std::vector<
            std::unique_ptr<MolecularOrbitalRenderer>
        > h2oOrbitalRenderers;

        for (
            const MolecularOrbitalIsosurface::Surface& surface :
            h2oSurfaces
        )
        {
            if (surface.empty())
            {
                h2oOrbitalRenderers.push_back(
                    nullptr
                );

                continue;
            }

            auto renderer =
                std::make_unique<MolecularOrbitalRenderer>();

            renderer->initialize();

            renderer->setSurface(
                surface
            );

            renderer->setPositiveColor(
                glm::vec3(
                    POSITIVE_RED,
                    POSITIVE_GREEN,
                    POSITIVE_BLUE
                )
            );

            renderer->setNegativeColor(
                glm::vec3(
                    NEGATIVE_RED,
                    NEGATIVE_GREEN,
                    NEGATIVE_BLUE
                )
            );

            renderer->setModelMatrix(
                glm::translate(
                    glm::mat4(1.0f),
                    glm::vec3(
                        H2O_POSITION_X,
                        H2O_POSITION_Y,
                        0.0f
                    )
                )
            );

            h2oOrbitalRenderers.push_back(
                std::move(renderer)
            );
        }


        // =====================================================
        // H2 NUCLEI
        // =====================================================

        MolecularNucleusRenderer h2NucleusRenderer;

        h2NucleusRenderer.initialize();

        h2NucleusRenderer.setMolecule(
            h2
        );

        h2NucleusRenderer.setColor(
            glm::vec3(
                NUCLEUS_RED,
                NUCLEUS_GREEN,
                NUCLEUS_BLUE
            )
        );

        h2NucleusRenderer.setModelMatrix(
            glm::translate(
                glm::mat4(1.0f),
                glm::vec3(
                    H2_POSITION_X,
                    H2_POSITION_Y,
                    0.0f
                )
            )
        );


        // =====================================================
        // H2O NUCLEI
        // =====================================================

        MolecularNucleusRenderer h2oNucleusRenderer;

        h2oNucleusRenderer.initialize();

        h2oNucleusRenderer.setMolecule(
            h2o
        );

        h2oNucleusRenderer.setColor(
            glm::vec3(
                NUCLEUS_RED,
                NUCLEUS_GREEN,
                NUCLEUS_BLUE
            )
        );

        h2oNucleusRenderer.setModelMatrix(
            glm::translate(
                glm::mat4(1.0f),
                glm::vec3(
                    H2O_POSITION_X,
                    H2O_POSITION_Y,
                    0.0f
                )
            )
        );


        // =====================================================
        // RENDER LOOP
        // =====================================================

        while (!navigation.shouldClose())
        {
            navigation.update();

            glClear(
                GL_COLOR_BUFFER_BIT |
                GL_DEPTH_BUFFER_BIT
            );

            const glm::mat4& view =
                navigation.getViewMatrix();

            const glm::mat4& projection =
                navigation.getProjectionMatrix();


            // -------------------------------------------------
            // Grid
            // -------------------------------------------------

            gridRenderer.render(
                view,
                projection
            );


            // -------------------------------------------------
            // H2 - ALL ORBITALS
            // -------------------------------------------------

            for (
                std::unique_ptr<MolecularOrbitalRenderer>& renderer :
                h2OrbitalRenderers
            )
            {
                if (!renderer)
                {
                    continue;
                }

                renderer->setViewMatrix(
                    view
                );

                renderer->setProjectionMatrix(
                    projection
                );

                renderer->render();
            }


            // -------------------------------------------------
            // H2 NUCLEI
            // -------------------------------------------------

            h2NucleusRenderer.setViewMatrix(
                view
            );

            h2NucleusRenderer.setProjectionMatrix(
                projection
            );

            h2NucleusRenderer.render();


            // -------------------------------------------------
            // H2O - ALL ORBITALS
            // -------------------------------------------------

            for (
                std::unique_ptr<MolecularOrbitalRenderer>& renderer :
                h2oOrbitalRenderers
            )
            {
                if (!renderer)
                {
                    continue;
                }

                renderer->setViewMatrix(
                    view
                );

                renderer->setProjectionMatrix(
                    projection
                );

                renderer->render();
            }


            // -------------------------------------------------
            // H2O NUCLEI
            // -------------------------------------------------

            h2oNucleusRenderer.setViewMatrix(
                view
            );

            h2oNucleusRenderer.setProjectionMatrix(
                projection
            );

            h2oNucleusRenderer.render();


            navigation.present();
        }
    }
    catch (const std::exception& exception)
    {
        std::cerr
            << "ERROR: "
            << exception.what()
            << '\n';

        return 1;
    }

    return 0;
}