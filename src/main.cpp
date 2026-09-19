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

#include <algorithm>
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
        maximum = std::max(
            maximum,
            std::abs(value)
        );
    }

    return maximum;
}


bool sameSpatialOrbital(
    const MolecularOrbital& first,
    const MolecularOrbital& second
)
{
    if (first.psi.size() != second.psi.size())
    {
        return false;
    }

    if (
        std::abs(
            first.eigenvalue -
            second.eigenvalue
        ) > 1.0e-10
    )
    {
        return false;
    }

    for (std::size_t i = 0; i < first.psi.size(); ++i)
    {
        if (
            std::abs(
                first.psi[i] -
                second.psi[i]
            ) > 1.0e-10
        )
        {
            return false;
        }
    }

    return true;
}


std::vector<MolecularOrbital> extractUniqueSpatialOrbitals(
    const std::vector<MolecularOrbital>& orbitals
)
{
    std::vector<MolecularOrbital> uniqueOrbitals;

    for (const MolecularOrbital& orbital : orbitals)
    {
        bool alreadyPresent = false;

        for (
            const MolecularOrbital& existing :
            uniqueOrbitals
        )
        {
            if (
                sameSpatialOrbital(
                    orbital,
                    existing
                )
            )
            {
                alreadyPresent = true;
                break;
            }
        }

        if (!alreadyPresent)
        {
            uniqueOrbitals.push_back(
                orbital
            );
        }
    }

    return uniqueOrbitals;
}


MolecularOrbitalIsosurface::Surface generateSurface(
    const CartesianGrid& grid,
    const MolecularOrbital& orbital
)
{
    return MolecularOrbitalIsosurface::generate(
        grid,
        orbital.psi,
        ISOVALUE
    );
}


std::unique_ptr<MolecularOrbitalRenderer>
createOrbitalRenderer(
    const MolecularOrbitalIsosurface::Surface& surface,
    float x,
    float y
)
{
    if (surface.empty())
    {
        return nullptr;
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
                x,
                y,
                0.0f
            )
        )
    );

    return renderer;
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
        // UNIQUE SPATIAL ORBITALS
        // =====================================================

        const std::vector<MolecularOrbital> h2Orbitals =
            extractUniqueSpatialOrbitals(
                h2Result.scf.molecularOrbitals
            );

        const std::vector<MolecularOrbital> h2oOrbitals =
            extractUniqueSpatialOrbitals(
                h2oResult.scf.molecularOrbitals
            );


        // =====================================================
        // INFORMATION
        // =====================================================

        std::cout
            << "\n========================================\n"
            << "H2 MOLECULAR ORBITALS\n"
            << "========================================\n";

        std::cout
            << "Returned orbitals: "
            << h2Result.scf.molecularOrbitals.size()
            << '\n';

        std::cout
            << "Unique spatial orbitals: "
            << h2Orbitals.size()
            << '\n';

        for (
            std::size_t i = 0;
            i < h2Orbitals.size();
            ++i
        )
        {
            const MolecularOrbital& orbital =
                h2Orbitals[i];

            std::cout
                << "MO "
                << i
                << " | eigenvalue = "
                << orbital.eigenvalue
                << " | electrons = "
                << orbital.electrons
                << " | max|psi| = "
                << maximumAbsoluteValue(
                    orbital.psi
                )
                << '\n';
        }


        std::cout
            << "\n========================================\n"
            << "H2O MOLECULAR ORBITALS\n"
            << "========================================\n";

        std::cout
            << "Returned orbitals: "
            << h2oResult.scf.molecularOrbitals.size()
            << '\n';

        std::cout
            << "Unique spatial orbitals: "
            << h2oOrbitals.size()
            << '\n';

        for (
            std::size_t i = 0;
            i < h2oOrbitals.size();
            ++i
        )
        {
            const MolecularOrbital& orbital =
                h2oOrbitals[i];

            std::cout
                << "MO "
                << i
                << " | eigenvalue = "
                << orbital.eigenvalue
                << " | electrons = "
                << orbital.electrons
                << " | max|psi| = "
                << maximumAbsoluteValue(
                    orbital.psi
                )
                << '\n';
        }


        // =====================================================
        // ISOSURFACES
        // =====================================================

        std::vector<
            MolecularOrbitalIsosurface::Surface
        > h2Surfaces;

        h2Surfaces.reserve(
            h2Orbitals.size()
        );

        for (
            const MolecularOrbital& orbital :
            h2Orbitals
        )
        {
            h2Surfaces.push_back(
                generateSurface(
                    h2Grid,
                    orbital
                )
            );
        }


        std::vector<
            MolecularOrbitalIsosurface::Surface
        > h2oSurfaces;

        h2oSurfaces.reserve(
            h2oOrbitals.size()
        );

        for (
            const MolecularOrbital& orbital :
            h2oOrbitals
        )
        {
            h2oSurfaces.push_back(
                generateSurface(
                    h2oGrid,
                    orbital
                )
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
        // H2 ORBITALS
        // =====================================================

        std::vector<
            std::unique_ptr<MolecularOrbitalRenderer>
        > h2OrbitalRenderers;

        for (
            const MolecularOrbitalIsosurface::Surface& surface :
            h2Surfaces
        )
        {
            h2OrbitalRenderers.push_back(
                createOrbitalRenderer(
                    surface,
                    H2_POSITION_X,
                    H2_POSITION_Y
                )
            );
        }


        // =====================================================
        // H2O ORBITALS
        // =====================================================

        std::vector<
            std::unique_ptr<MolecularOrbitalRenderer>
        > h2oOrbitalRenderers;

        for (
            const MolecularOrbitalIsosurface::Surface& surface :
            h2oSurfaces
        )
        {
            h2oOrbitalRenderers.push_back(
                createOrbitalRenderer(
                    surface,
                    H2O_POSITION_X,
                    H2O_POSITION_Y
                )
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
            // GRID
            // -------------------------------------------------

            gridRenderer.render(
                view,
                projection
            );


            // -------------------------------------------------
            // H2
            // -------------------------------------------------

            for (
                std::unique_ptr<MolecularOrbitalRenderer>&
                    renderer :
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


            h2NucleusRenderer.setViewMatrix(
                view
            );

            h2NucleusRenderer.setProjectionMatrix(
                projection
            );

            h2NucleusRenderer.render();


            // -------------------------------------------------
            // H2O
            // -------------------------------------------------

            for (
                std::unique_ptr<MolecularOrbitalRenderer>&
                    renderer :
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