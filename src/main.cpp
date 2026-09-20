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
#include <string>
#include <vector>

namespace
{

constexpr int WINDOW_WIDTH = 1280;
constexpr int WINDOW_HEIGHT = 720;

constexpr int GRID_POINTS = 9;

constexpr double GRID_MIN = -3.8;
constexpr double GRID_MAX = 4.2;

constexpr double ISOVALUE = 0.04;


// =============================================================
// ORBITAL POSITIONS
// =============================================================

constexpr float N2_X = -6.0f;
constexpr float N2_Y = 2.5f;

constexpr float O2_X = 0.0f;
constexpr float O2_Y = 2.5f;

constexpr float H2O_X = 6.0f;
constexpr float H2O_Y = 2.5f;

constexpr float N2_SMOOTH_X = -6.0f;
constexpr float N2_SMOOTH_Y = -3.5f;

constexpr float O2_SMOOTH_X = 0.0f;
constexpr float O2_SMOOTH_Y = -3.5f;

constexpr float CO2_X = 6.0f;
constexpr float CO2_Y = -3.5f;


// =============================================================
// DENSITY POSITIONS
// =============================================================

constexpr float N2_DENSITY_X = -6.0f;
constexpr float N2_DENSITY_Y = -8.0f;

constexpr float O2_DENSITY_X = 0.0f;
constexpr float O2_DENSITY_Y = -8.0f;

constexpr float H2O_DENSITY_X = 6.0f;
constexpr float H2O_DENSITY_Y = -8.0f;

constexpr float CO2_DENSITY_X = 6.0f;
constexpr float CO2_DENSITY_Y = -8.0f;


// =============================================================
// COLORS
// =============================================================

constexpr float POSITIVE_RED = 0.15f;
constexpr float POSITIVE_GREEN = 0.45f;
constexpr float POSITIVE_BLUE = 1.0f;

constexpr float NEGATIVE_RED = 1.0f;
constexpr float NEGATIVE_GREEN = 0.20f;
constexpr float NEGATIVE_BLUE = 0.20f;

constexpr float DENSITY_RED = 0.20f;
constexpr float DENSITY_GREEN = 0.75f;
constexpr float DENSITY_BLUE = 0.30f;

constexpr float NUCLEUS_RED = 0.90f;
constexpr float NUCLEUS_GREEN = 0.90f;
constexpr float NUCLEUS_BLUE = 0.90f;


// =============================================================
// DATA
// =============================================================

struct MoleculeData
{
    std::string name;
    Molecule molecule;
    CartesianGrid grid;
    int multiplicity;

    MolecularResult result;

    std::vector<MolecularOrbital> orbitals;
};


// =============================================================
// GRID
// =============================================================

CartesianGrid createGrid()
{
    return CartesianGrid(
        GRID_POINTS,
        GRID_POINTS,
        GRID_POINTS,
        GRID_MIN,
        GRID_MAX,
        GRID_MIN,
        GRID_MAX,
        GRID_MIN,
        GRID_MAX
    );
}


// =============================================================
// WAVEFUNCTION
// =============================================================

std::vector<double> squareWavefunction(
    const std::vector<double>& psi
)
{
    std::vector<double> density(
        psi.size()
    );

    for (std::size_t i = 0; i < psi.size(); ++i)
    {
        density[i] =
            psi[i] * psi[i];
    }

    return density;
}


// =============================================================
// MATHEMATICAL SMOOTHING
// =============================================================

std::vector<double> smoothWavefunction(
    const CartesianGrid& grid,
    const std::vector<double>& psi
)
{
    if (psi.size() != grid.getSize())
    {
        return psi;
    }

    std::vector<double> smoothed(
        psi.size(),
        0.0
    );

    const std::size_t nx =
        grid.getNx();

    const std::size_t ny =
        grid.getNy();

    const std::size_t nz =
        grid.getNz();

    for (std::size_t k = 0; k < nz; ++k)
    {
        for (std::size_t j = 0; j < ny; ++j)
        {
            for (std::size_t i = 0; i < nx; ++i)
            {
                double sum = 0.0;
                std::size_t count = 0;

                for (int dk = -1; dk <= 1; ++dk)
                {
                    for (int dj = -1; dj <= 1; ++dj)
                    {
                        for (int di = -1; di <= 1; ++di)
                        {
                            const int ii =
                                static_cast<int>(i) + di;

                            const int jj =
                                static_cast<int>(j) + dj;

                            const int kk =
                                static_cast<int>(k) + dk;

                            if (ii < 0 ||
                                jj < 0 ||
                                kk < 0 ||
                                ii >= static_cast<int>(nx) ||
                                jj >= static_cast<int>(ny) ||
                                kk >= static_cast<int>(nz))
                            {
                                continue;
                            }

                            const std::size_t index =
                                grid.getIndex(
                                    static_cast<std::size_t>(ii),
                                    static_cast<std::size_t>(jj),
                                    static_cast<std::size_t>(kk)
                                );

                            sum += psi[index];

                            ++count;
                        }
                    }
                }

                const std::size_t index =
                    grid.getIndex(i, j, k);

                if (count > 0)
                {
                    smoothed[index] =
                        sum /
                        static_cast<double>(count);
                }
                else
                {
                    smoothed[index] =
                        psi[index];
                }
            }
        }
    }

    return smoothed;
}


// =============================================================
// UNIQUE SPATIAL ORBITALS
// =============================================================

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


// =============================================================
// SURFACES
// =============================================================

MolecularOrbitalIsosurface::Surface generateOrbitalSurface(
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


MolecularOrbitalIsosurface::Surface generateSmoothOrbitalSurface(
    const CartesianGrid& grid,
    const MolecularOrbital& orbital
)
{
    const std::vector<double> smoothedPsi =
        smoothWavefunction(
            grid,
            orbital.psi
        );

    return MolecularOrbitalIsosurface::generate(
        grid,
        smoothedPsi,
        ISOVALUE
    );
}


MolecularOrbitalIsosurface::Surface generateDensitySurface(
    const CartesianGrid& grid,
    const MolecularOrbital& orbital
)
{
    const std::vector<double> density =
        squareWavefunction(
            orbital.psi
        );

    return MolecularOrbitalIsosurface::generate(
        grid,
        density,
        ISOVALUE
    );
}


// =============================================================
// ORBITAL RENDERER
// =============================================================

std::unique_ptr<MolecularOrbitalRenderer> createOrbitalRenderer(
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


// =============================================================
// DENSITY RENDERER
// =============================================================

std::unique_ptr<MolecularOrbitalRenderer> createDensityRenderer(
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
            DENSITY_RED,
            DENSITY_GREEN,
            DENSITY_BLUE
        )
    );

    renderer->setNegativeColor(
        glm::vec3(
            DENSITY_RED,
            DENSITY_GREEN,
            DENSITY_BLUE
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


// =============================================================
// RENDERER LIST
// =============================================================

using RendererList =
    std::vector<
        std::unique_ptr<MolecularOrbitalRenderer>
    >;


void renderOrbitalList(
    RendererList& renderers,
    const glm::mat4& view,
    const glm::mat4& projection
)
{
    for (
        std::unique_ptr<MolecularOrbitalRenderer>& renderer :
        renderers
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
}


// =============================================================
// NUCLEUS RENDERER
// =============================================================

void configureNucleusRenderer(
    MolecularNucleusRenderer& renderer,
    const Molecule& molecule,
    float x,
    float y
)
{
    renderer.initialize();

    renderer.setMolecule(
        molecule
    );

    renderer.setColor(
        glm::vec3(
            NUCLEUS_RED,
            NUCLEUS_GREEN,
            NUCLEUS_BLUE
        )
    );

    renderer.setModelMatrix(
        glm::translate(
            glm::mat4(1.0f),
            glm::vec3(
                x,
                y,
                0.0f
            )
        )
    );
}


void renderNucleus(
    MolecularNucleusRenderer& renderer,
    const glm::mat4& view,
    const glm::mat4& projection
)
{
    renderer.setViewMatrix(
        view
    );

    renderer.setProjectionMatrix(
        projection
    );

    renderer.render();
}

}


// =============================================================
// MAIN
// =============================================================

int main()
{
    try
    {
        PZ81 functional;


        // =====================================================
        // N2
        // =====================================================

        Molecule n2(0);

        n2.addNucleus(
            7,
            -1.037,
            0.0,
            0.0
        );

        n2.addNucleus(
            7,
            1.037,
            0.0,
            0.0
        );


        // =====================================================
        // O2
        // MULTIPLICITY = 3
        // =====================================================

        Molecule o2(0);

        o2.addNucleus(
            8,
            -1.14,
            0.0,
            0.0
        );

        o2.addNucleus(
            8,
            1.14,
            0.0,
            0.0
        );


        // =====================================================
        // H2O
        // =====================================================

        Molecule h2o(0);

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


        // =====================================================
        // CO2
        // =====================================================

        Molecule co2(0);

        co2.addNucleus(
            6,
            0.0,
            0.0,
            0.0
        );

        co2.addNucleus(
            8,
            -2.192,
            0.0,
            0.0
        );

        co2.addNucleus(
            8,
            2.192,
            0.0,
            0.0
        );


        // =====================================================
        // DATA
        // =====================================================

        MoleculeData n2Data{
            "N2",
            n2,
            createGrid(),
            1,
            {},
            {}
        };

        MoleculeData o2Data{
            "O2",
            o2,
            createGrid(),
            3,
            {},
            {}
        };

        MoleculeData h2oData{
            "H2O",
            h2o,
            createGrid(),
            1,
            {},
            {}
        };

        MoleculeData co2Data{
            "CO2",
            co2,
            createGrid(),
            1,
            {},
            {}
        };


        std::vector<MoleculeData*> molecules{
            &n2Data,
            &o2Data,
            &h2oData,
            &co2Data
        };


        // =====================================================
        // DFT
        // =====================================================

        for (MoleculeData* data : molecules)
        {
            std::cout
                << "\nCalculating "
                << data->name
                << "...\n";

            data->result =
                solveMolecularSelfConsistentField(
                    data->grid,
                    data->molecule,
                    data->molecule.getCharge(),
                    data->multiplicity,
                    functional
                );

            std::cout
                << data->name
                << " | Converged = "
                << (
                    data->result.scf.converged
                        ? "Yes"
                        : "No"
                )
                << " | Iterations = "
                << data->result.scf.iterations
                << " | E = "
                << data->result.scf.energy.total
                << " Ha\n";

            if (!data->result.scf.converged)
            {
                continue;
            }

            data->orbitals =
                extractUniqueSpatialOrbitals(
                    data->result.scf.molecularOrbitals
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
        // ORBITAL RENDERERS
        // =====================================================

        RendererList n2OrbitalRenderers;
        RendererList o2OrbitalRenderers;
        RendererList h2oOrbitalRenderers;
        RendererList co2OrbitalRenderers;

        RendererList n2SmoothOrbitalRenderers;
        RendererList o2SmoothOrbitalRenderers;

        RendererList n2DensityRenderers;
        RendererList o2DensityRenderers;
        RendererList h2oDensityRenderers;
        RendererList co2DensityRenderers;


        // =====================================================
        // N2
        // ORIGINAL + SMOOTH
        // =====================================================

        for (const MolecularOrbital& orbital : n2Data.orbitals)
        {
            n2OrbitalRenderers.push_back(
                createOrbitalRenderer(
                    generateOrbitalSurface(
                        n2Data.grid,
                        orbital
                    ),
                    N2_X,
                    N2_Y
                )
            );

            n2SmoothOrbitalRenderers.push_back(
                createOrbitalRenderer(
                    generateSmoothOrbitalSurface(
                        n2Data.grid,
                        orbital
                    ),
                    N2_SMOOTH_X,
                    N2_SMOOTH_Y
                )
            );

            n2DensityRenderers.push_back(
                createDensityRenderer(
                    generateDensitySurface(
                        n2Data.grid,
                        orbital
                    ),
                    N2_DENSITY_X,
                    N2_DENSITY_Y
                )
            );
        }


        // =====================================================
        // O2
        // ORIGINAL + SMOOTH
        // =====================================================

        for (const MolecularOrbital& orbital : o2Data.orbitals)
        {
            o2OrbitalRenderers.push_back(
                createOrbitalRenderer(
                    generateOrbitalSurface(
                        o2Data.grid,
                        orbital
                    ),
                    O2_X,
                    O2_Y
                )
            );

            o2SmoothOrbitalRenderers.push_back(
                createOrbitalRenderer(
                    generateSmoothOrbitalSurface(
                        o2Data.grid,
                        orbital
                    ),
                    O2_SMOOTH_X,
                    O2_SMOOTH_Y
                )
            );

            o2DensityRenderers.push_back(
                createDensityRenderer(
                    generateDensitySurface(
                        o2Data.grid,
                        orbital
                    ),
                    O2_DENSITY_X,
                    O2_DENSITY_Y
                )
            );
        }


        // =====================================================
        // H2O
        // =====================================================

        for (const MolecularOrbital& orbital : h2oData.orbitals)
        {
            h2oOrbitalRenderers.push_back(
                createOrbitalRenderer(
                    generateOrbitalSurface(
                        h2oData.grid,
                        orbital
                    ),
                    H2O_X,
                    H2O_Y
                )
            );

            h2oDensityRenderers.push_back(
                createDensityRenderer(
                    generateDensitySurface(
                        h2oData.grid,
                        orbital
                    ),
                    H2O_DENSITY_X,
                    H2O_DENSITY_Y
                )
            );
        }


        // =====================================================
        // CO2
        // =====================================================

        for (const MolecularOrbital& orbital : co2Data.orbitals)
        {
            co2OrbitalRenderers.push_back(
                createOrbitalRenderer(
                    generateOrbitalSurface(
                        co2Data.grid,
                        orbital
                    ),
                    CO2_X,
                    CO2_Y
                )
            );

            co2DensityRenderers.push_back(
                createDensityRenderer(
                    generateDensitySurface(
                        co2Data.grid,
                        orbital
                    ),
                    CO2_DENSITY_X,
                    CO2_DENSITY_Y
                )
            );
        }


        // =====================================================
        // NUCLEI
        // =====================================================

        MolecularNucleusRenderer n2NucleusRenderer;
        MolecularNucleusRenderer o2NucleusRenderer;
        MolecularNucleusRenderer h2oNucleusRenderer;
        MolecularNucleusRenderer co2NucleusRenderer;

        configureNucleusRenderer(
            n2NucleusRenderer,
            n2,
            N2_X,
            N2_Y
        );

        configureNucleusRenderer(
            o2NucleusRenderer,
            o2,
            O2_X,
            O2_Y
        );

        configureNucleusRenderer(
            h2oNucleusRenderer,
            h2o,
            H2O_X,
            H2O_Y
        );

        configureNucleusRenderer(
            co2NucleusRenderer,
            co2,
            CO2_X,
            CO2_Y
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
            // N2 ORIGINAL
            // -------------------------------------------------

            renderOrbitalList(
                n2OrbitalRenderers,
                view,
                projection
            );

            renderNucleus(
                n2NucleusRenderer,
                view,
                projection
            );


            // -------------------------------------------------
            // O2 ORIGINAL
            // -------------------------------------------------

            renderOrbitalList(
                o2OrbitalRenderers,
                view,
                projection
            );

            renderNucleus(
                o2NucleusRenderer,
                view,
                projection
            );


            // -------------------------------------------------
            // H2O
            // -------------------------------------------------

            renderOrbitalList(
                h2oOrbitalRenderers,
                view,
                projection
            );

            renderNucleus(
                h2oNucleusRenderer,
                view,
                projection
            );


            // -------------------------------------------------
            // N2 SMOOTH
            // -------------------------------------------------

            renderOrbitalList(
                n2SmoothOrbitalRenderers,
                view,
                projection
            );


            // -------------------------------------------------
            // O2 SMOOTH
            // -------------------------------------------------

            renderOrbitalList(
                o2SmoothOrbitalRenderers,
                view,
                projection
            );


            // -------------------------------------------------
            // CO2
            // -------------------------------------------------

            renderOrbitalList(
                co2OrbitalRenderers,
                view,
                projection
            );

            renderNucleus(
                co2NucleusRenderer,
                view,
                projection
            );


            // -------------------------------------------------
            // DENSITIES
            // -------------------------------------------------

            renderOrbitalList(
                n2DensityRenderers,
                view,
                projection
            );

            renderOrbitalList(
                o2DensityRenderers,
                view,
                projection
            );

            renderOrbitalList(
                h2oDensityRenderers,
                view,
                projection
            );

            renderOrbitalList(
                co2DensityRenderers,
                view,
                projection
            );


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