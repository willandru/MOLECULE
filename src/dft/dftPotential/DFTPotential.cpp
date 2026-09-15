#include "DFTPotential.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace
{
    constexpr double PI =
        3.1415926535897932384626433832795;

    constexpr double DEFAULT_SOFTENING =
        0.15;

    constexpr double MIN_DENSITY =
        1.0e-14;

    constexpr double MIN_DISTANCE =
        1.0e-12;

    constexpr double POISSON_TOLERANCE =
        1.0e-8;

    constexpr int PRE_SMOOTH =
        2;

    constexpr int POST_SMOOTH =
        2;

    constexpr int COARSE_SMOOTH =
        64;
}


// ================================================================
// CONSTRUCTORES
// ================================================================

DFTPotential::DFTPotential()
    :
    grid(nullptr),
    softening(DEFAULT_SOFTENING),
    poissonResidual(
        std::numeric_limits<double>::infinity()
    ),
    poissonCycles(0),
    poissonConverged(false)
{
}


DFTPotential::DFTPotential(
    const DFTGrid& grid
)
    :
    grid(nullptr),
    softening(DEFAULT_SOFTENING),
    poissonResidual(
        std::numeric_limits<double>::infinity()
    ),
    poissonCycles(0),
    poissonConverged(false)
{
    initialize(grid);
}


// ================================================================
// INICIALIZACIÓN
// ================================================================

void DFTPotential::initialize(
    const DFTGrid& grid
)
{
    if (grid.getPointCount() == 0)
    {
        throw std::invalid_argument(
            "DFTPotential: grid has no points."
        );
    }

    if (grid.getNx() < 3 ||
        grid.getNy() < 3 ||
        grid.getNz() < 3)
    {
        throw std::invalid_argument(
            "DFTPotential: grid dimensions must be at least 3 in every direction."
        );
    }

    if (!std::isfinite(grid.getSpacing()) ||
        grid.getSpacing() <= 0.0)
    {
        throw std::invalid_argument(
            "DFTPotential: grid spacing must be positive and finite."
        );
    }

    this->grid = &grid;

    const std::size_t count =
        grid.getPointCount();

    externalPotential.assign(
        count,
        0.0
    );

    hartreePotential.assign(
        count,
        0.0
    );

    exchangePotential.assign(
        count,
        0.0
    );

    poissonResidual =
        std::numeric_limits<double>::infinity();

    poissonCycles =
        0;

    poissonConverged =
        false;
}


// ================================================================
// POTENCIAL EXTERNO
// ================================================================

void DFTPotential::calculateExternalPotential(
    const std::vector<int>& nuclearCharges,
    const std::vector<glm::dvec3>& nuclearPositions
)
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "DFTPotential: not initialized."
        );
    }

    if (nuclearCharges.size() !=
        nuclearPositions.size())
    {
        throw std::invalid_argument(
            "DFTPotential: nuclear charges and positions have different sizes."
        );
    }

    std::fill(
        externalPotential.begin(),
        externalPotential.end(),
        0.0
    );

    const std::size_t count =
        grid->getPointCount();

    for (std::size_t i = 0;
         i < count;
         ++i)
    {
        const glm::dvec3 r =
            grid->getPosition(i);

        double value =
            0.0;

        for (std::size_t nucleus = 0;
             nucleus < nuclearCharges.size();
             ++nucleus)
        {
            const int Z =
                nuclearCharges[nucleus];

            if (Z <= 0)
            {
                throw std::invalid_argument(
                    "DFTPotential: nuclear charge must be positive."
                );
            }

            const glm::dvec3 difference =
                r - nuclearPositions[nucleus];

            const double r2 =
                glm::dot(
                    difference,
                    difference
                );

            const double denominator =
                std::sqrt(
                    r2 +
                    softening * softening
                );

            value -=
                static_cast<double>(Z) /
                denominator;
        }

        externalPotential[i] =
            value;
    }
}


// ================================================================
// HARTREE
// ================================================================

void DFTPotential::calculateHartreePotential(
    const DFTDensity& density,
    int iterations
)
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "DFTPotential: not initialized."
        );
    }

    if (iterations <= 0)
    {
        throw std::invalid_argument(
            "DFTPotential: Hartree maximum cycles must be positive."
        );
    }

    if (density.getGrid().getPointCount() !=
        grid->getPointCount())
    {
        throw std::invalid_argument(
            "DFTPotential: incompatible density grid."
        );
    }

    solvePoisson(
        density,
        hartreePotential,
        iterations
    );
}


// ================================================================
// LAPLACIAN
// ================================================================

double DFTPotential::calculateLaplacian(
    const std::vector<double>& potential,
    int x,
    int y,
    int z
) const
{
    const int nx =
        grid->getNx();

    const int ny =
        grid->getNy();

    const int nz =
        grid->getNz();

    if (x <= 0 ||
        x >= nx - 1 ||
        y <= 0 ||
        y >= ny - 1 ||
        z <= 0 ||
        z >= nz - 1)
    {
        return 0.0;
    }

    const double h =
        grid->getSpacing();

    const double h2 =
        h * h;

    const std::size_t center =
        grid->getIndex(
            x,
            y,
            z
        );

    const double neighbors =
        potential[
            grid->getIndex(
                x + 1,
                y,
                z
            )
        ]
        +
        potential[
            grid->getIndex(
                x - 1,
                y,
                z
            )
        ]
        +
        potential[
            grid->getIndex(
                x,
                y + 1,
                z
            )
        ]
        +
        potential[
            grid->getIndex(
                x,
                y - 1,
                z
            )
        ]
        +
        potential[
            grid->getIndex(
                x,
                y,
                z + 1
            )
        ]
        +
        potential[
            grid->getIndex(
                x,
                y,
                z - 1
            )
        ];

    return
        (
            neighbors -
            6.0 * potential[center]
        )
        /
        h2;
}


// ================================================================
// CARGA ELECTRÓNICA TOTAL
// ================================================================

double DFTPotential::calculateDensityCharge(
    const DFTDensity& density
) const
{
    return density.calculateElectronCount();
}


// ================================================================
// CONDICIONES DE FRONTERA
// ================================================================

void DFTPotential::applyBoundaryConditions(
    std::vector<double>& potential,
    const DFTDensity& density
) const
{
    const int nx =
        grid->getNx();

    const int ny =
        grid->getNy();

    const int nz =
        grid->getNz();

    const double electronCount =
        calculateDensityCharge(
            density
        );

    if (!std::isfinite(electronCount) ||
        electronCount <= 0.0)
    {
        return;
    }

    const double dV =
        grid->getVolumeElement();

    glm::dvec3 center(
        0.0
    );

    double totalCharge =
        0.0;

    for (std::size_t i = 0;
         i < grid->getPointCount();
         ++i)
    {
        const double rho =
            density.get(i);

        if (!std::isfinite(rho))
        {
            throw std::runtime_error(
                "DFTPotential: density contains a non-finite value."
            );
        }

        center +=
            rho *
            grid->getPosition(i) *
            dV;

        totalCharge +=
            rho *
            dV;
    }

    if (totalCharge > MIN_DENSITY)
    {
        center /=
            totalCharge;
    }

    auto setBoundary =
        [&](int x, int y, int z)
        {
            const std::size_t index =
                grid->getIndex(
                    x,
                    y,
                    z
                );

            const glm::dvec3 difference =
                grid->getPosition(index) -
                center;

            const double distance =
                glm::length(
                    difference
                );

            if (distance > MIN_DISTANCE)
            {
                potential[index] =
                    electronCount /
                    distance;
            }
        };

    for (int x = 0; x < nx; ++x)
    {
        for (int y = 0; y < ny; ++y)
        {
            setBoundary(
                x,
                y,
                0
            );

            setBoundary(
                x,
                y,
                nz - 1
            );
        }
    }

    for (int x = 0; x < nx; ++x)
    {
        for (int z = 0; z < nz; ++z)
        {
            setBoundary(
                x,
                0,
                z
            );

            setBoundary(
                x,
                ny - 1,
                z
            );
        }
    }

    for (int y = 0; y < ny; ++y)
    {
        for (int z = 0; z < nz; ++z)
        {
            setBoundary(
                0,
                y,
                z
            );

            setBoundary(
                nx - 1,
                y,
                z
            );
        }
    }
}


// ================================================================
// SUAVIZADO GAUSS-SEIDEL
// ================================================================

void DFTPotential::smooth(
    std::vector<double>& solution,
    const std::vector<double>& rhs,
    int nx,
    int ny,
    int nz,
    double hx,
    double hy,
    double hz,
    int iterations
) const
{
    if (nx < 3 ||
        ny < 3 ||
        nz < 3)
    {
        return;
    }

    const double inverseHx2 =
        1.0 / (hx * hx);

    const double inverseHy2 =
        1.0 / (hy * hy);

    const double inverseHz2 =
        1.0 / (hz * hz);

    const double diagonal =
        2.0 * (
            inverseHx2 +
            inverseHy2 +
            inverseHz2
        );

    for (int iteration = 0;
         iteration < iterations;
         ++iteration)
    {
        for (int z = 1;
             z < nz - 1;
             ++z)
        {
            for (int y = 1;
                 y < ny - 1;
                 ++y)
            {
                for (int x = 1;
                     x < nx - 1;
                     ++x)
                {
                    const std::size_t index =
                        static_cast<std::size_t>(x)
                        +
                        static_cast<std::size_t>(nx)
                        *
                        (
                            static_cast<std::size_t>(y)
                            +
                            static_cast<std::size_t>(ny)
                            *
                            static_cast<std::size_t>(z)
                        );

                    const std::size_t xp =
                        static_cast<std::size_t>(x + 1)
                        +
                        static_cast<std::size_t>(nx)
                        *
                        (
                            static_cast<std::size_t>(y)
                            +
                            static_cast<std::size_t>(ny)
                            *
                            static_cast<std::size_t>(z)
                        );

                    const std::size_t xm =
                        static_cast<std::size_t>(x - 1)
                        +
                        static_cast<std::size_t>(nx)
                        *
                        (
                            static_cast<std::size_t>(y)
                            +
                            static_cast<std::size_t>(ny)
                            *
                            static_cast<std::size_t>(z)
                        );

                    const std::size_t yp =
                        static_cast<std::size_t>(x)
                        +
                        static_cast<std::size_t>(nx)
                        *
                        (
                            static_cast<std::size_t>(y + 1)
                            +
                            static_cast<std::size_t>(ny)
                            *
                            static_cast<std::size_t>(z)
                        );

                    const std::size_t ym =
                        static_cast<std::size_t>(x)
                        +
                        static_cast<std::size_t>(nx)
                        *
                        (
                            static_cast<std::size_t>(y - 1)
                            +
                            static_cast<std::size_t>(ny)
                            *
                            static_cast<std::size_t>(z)
                        );

                    const std::size_t zp =
                        static_cast<std::size_t>(x)
                        +
                        static_cast<std::size_t>(nx)
                        *
                        (
                            static_cast<std::size_t>(y)
                            +
                            static_cast<std::size_t>(ny)
                            *
                            static_cast<std::size_t>(z + 1)
                        );

                    const std::size_t zm =
                        static_cast<std::size_t>(x)
                        +
                        static_cast<std::size_t>(nx)
                        *
                        (
                            static_cast<std::size_t>(y)
                            +
                            static_cast<std::size_t>(ny)
                            *
                            static_cast<std::size_t>(z - 1)
                        );

                    const double neighborContribution =
                        (
                            solution[xp] +
                            solution[xm]
                        )
                        *
                        inverseHx2
                        +
                        (
                            solution[yp] +
                            solution[ym]
                        )
                        *
                        inverseHy2
                        +
                        (
                            solution[zp] +
                            solution[zm]
                        )
                        *
                        inverseHz2;

                    solution[index] =
                        (
                            neighborContribution -
                            rhs[index]
                        )
                        /
                        diagonal;
                }
            }
        }
    }
}


// ================================================================
// RESIDUO
// ================================================================

std::vector<double>
DFTPotential::calculateResidual(
    const std::vector<double>& solution,
    const std::vector<double>& rhs,
    int nx,
    int ny,
    int nz,
    double hx,
    double hy,
    double hz
) const
{
    const std::size_t count =
        static_cast<std::size_t>(nx)
        *
        static_cast<std::size_t>(ny)
        *
        static_cast<std::size_t>(nz);

    if (solution.size() != count ||
        rhs.size() != count)
    {
        throw std::invalid_argument(
            "DFTPotential: invalid Poisson vector size."
        );
    }

    std::vector<double> residual(
        count,
        0.0
    );

    const double inverseHx2 =
        1.0 / (hx * hx);

    const double inverseHy2 =
        1.0 / (hy * hy);

    const double inverseHz2 =
        1.0 / (hz * hz);

    for (int z = 1;
         z < nz - 1;
         ++z)
    {
        for (int y = 1;
             y < ny - 1;
             ++y)
        {
            for (int x = 1;
                 x < nx - 1;
                 ++x)
            {
                const std::size_t index =
                    static_cast<std::size_t>(x)
                    +
                    static_cast<std::size_t>(nx)
                    *
                    (
                        static_cast<std::size_t>(y)
                        +
                        static_cast<std::size_t>(ny)
                        *
                        static_cast<std::size_t>(z)
                    );

                const std::size_t xp =
                    static_cast<std::size_t>(x + 1)
                    +
                    static_cast<std::size_t>(nx)
                    *
                    (
                        static_cast<std::size_t>(y)
                        +
                        static_cast<std::size_t>(ny)
                        *
                        static_cast<std::size_t>(z)
                    );

                const std::size_t xm =
                    static_cast<std::size_t>(x - 1)
                    +
                    static_cast<std::size_t>(nx)
                    *
                    (
                        static_cast<std::size_t>(y)
                        +
                        static_cast<std::size_t>(ny)
                        *
                        static_cast<std::size_t>(z)
                    );

                const std::size_t yp =
                    static_cast<std::size_t>(x)
                    +
                    static_cast<std::size_t>(nx)
                    *
                    (
                        static_cast<std::size_t>(y + 1)
                        +
                        static_cast<std::size_t>(ny)
                        *
                        static_cast<std::size_t>(z)
                    );

                const std::size_t ym =
                    static_cast<std::size_t>(x)
                    +
                    static_cast<std::size_t>(nx)
                    *
                    (
                        static_cast<std::size_t>(y - 1)
                        +
                        static_cast<std::size_t>(ny)
                        *
                        static_cast<std::size_t>(z)
                    );

                const std::size_t zp =
                    static_cast<std::size_t>(x)
                    +
                    static_cast<std::size_t>(nx)
                    *
                    (
                        static_cast<std::size_t>(y)
                        +
                        static_cast<std::size_t>(ny)
                        *
                        static_cast<std::size_t>(z + 1)
                    );

                const std::size_t zm =
                    static_cast<std::size_t>(x)
                    +
                    static_cast<std::size_t>(nx)
                    *
                    (
                        static_cast<std::size_t>(y)
                        +
                        static_cast<std::size_t>(ny)
                        *
                        static_cast<std::size_t>(z - 1)
                    );

                const double laplacian =
                    (
                        solution[xp] -
                        2.0 * solution[index] +
                        solution[xm]
                    )
                    *
                    inverseHx2
                    +
                    (
                        solution[yp] -
                        2.0 * solution[index] +
                        solution[ym]
                    )
                    *
                    inverseHy2
                    +
                    (
                        solution[zp] -
                        2.0 * solution[index] +
                        solution[zm]
                    )
                    *
                    inverseHz2;

                residual[index] =
                    rhs[index] -
                    laplacian;
            }
        }
    }

    return residual;
}


// ================================================================
// NORMA DEL RESIDUO
// ================================================================

double DFTPotential::calculateResidualNorm(
    const std::vector<double>& residual,
    int nx,
    int ny,
    int nz,
    double hx,
    double hy,
    double hz
) const
{
    const std::size_t expectedSize =
        static_cast<std::size_t>(nx)
        *
        static_cast<std::size_t>(ny)
        *
        static_cast<std::size_t>(nz);

    if (residual.size() != expectedSize)
    {
        throw std::invalid_argument(
            "DFTPotential: invalid residual size."
        );
    }

    const double volumeElement =
        hx *
        hy *
        hz;

    double sum =
        0.0;

    for (int z = 1;
         z < nz - 1;
         ++z)
    {
        for (int y = 1;
             y < ny - 1;
             ++y)
        {
            for (int x = 1;
                 x < nx - 1;
                 ++x)
            {
                const std::size_t index =
                    static_cast<std::size_t>(x)
                    +
                    static_cast<std::size_t>(nx)
                    *
                    (
                        static_cast<std::size_t>(y)
                        +
                        static_cast<std::size_t>(ny)
                        *
                        static_cast<std::size_t>(z)
                    );

                const double value =
                    residual[index];

                sum +=
                    value *
                    value *
                    volumeElement;
            }
        }
    }

    return std::sqrt(
        std::max(
            sum,
            0.0
        )
    );
}


// ================================================================
// RESTRICTION
// ================================================================

std::vector<double>
DFTPotential::restrictResidual(
    const std::vector<double>& fineResidual,
    int fineNx,
    int fineNy,
    int fineNz,
    int coarseNx,
    int coarseNy,
    int coarseNz
) const
{
    const std::size_t fineCount =
        static_cast<std::size_t>(fineNx)
        *
        static_cast<std::size_t>(fineNy)
        *
        static_cast<std::size_t>(fineNz);

    if (fineResidual.size() != fineCount)
    {
        throw std::invalid_argument(
            "DFTPotential: invalid fine residual size."
        );
    }

    std::vector<double> coarseResidual(
        static_cast<std::size_t>(coarseNx)
        *
        static_cast<std::size_t>(coarseNy)
        *
        static_cast<std::size_t>(coarseNz),
        0.0
    );

    for (int z = 1;
         z < coarseNz - 1;
         ++z)
    {
        const double fineZ =
            static_cast<double>(z)
            *
            static_cast<double>(fineNz - 1)
            /
            static_cast<double>(coarseNz - 1);

        for (int y = 1;
             y < coarseNy - 1;
             ++y)
        {
            const double fineY =
                static_cast<double>(y)
                *
                static_cast<double>(fineNy - 1)
                /
                static_cast<double>(coarseNy - 1);

            for (int x = 1;
                 x < coarseNx - 1;
                 ++x)
            {
                const double fineX =
                    static_cast<double>(x)
                    *
                    static_cast<double>(fineNx - 1)
                    /
                    static_cast<double>(coarseNx - 1);

                const std::size_t coarseIndex =
                    static_cast<std::size_t>(x)
                    +
                    static_cast<std::size_t>(coarseNx)
                    *
                    (
                        static_cast<std::size_t>(y)
                        +
                        static_cast<std::size_t>(coarseNy)
                        *
                        static_cast<std::size_t>(z)
                    );

                coarseResidual[coarseIndex] =
                    sampleTrilinear(
                        fineResidual,
                        fineNx,
                        fineNy,
                        fineNz,
                        fineX,
                        fineY,
                        fineZ
                    );
            }
        }
    }

    return coarseResidual;
}


// ================================================================
// PROLONGATION
// ================================================================

void DFTPotential::prolongateAndAdd(
    const std::vector<double>& coarseCorrection,
    int coarseNx,
    int coarseNy,
    int coarseNz,
    std::vector<double>& fineSolution,
    int fineNx,
    int fineNy,
    int fineNz
) const
{
    const std::size_t coarseCount =
        static_cast<std::size_t>(coarseNx)
        *
        static_cast<std::size_t>(coarseNy)
        *
        static_cast<std::size_t>(coarseNz);

    const std::size_t fineCount =
        static_cast<std::size_t>(fineNx)
        *
        static_cast<std::size_t>(fineNy)
        *
        static_cast<std::size_t>(fineNz);

    if (coarseCorrection.size() != coarseCount ||
        fineSolution.size() != fineCount)
    {
        throw std::invalid_argument(
            "DFTPotential: invalid multigrid vector size."
        );
    }

    for (int z = 1;
         z < fineNz - 1;
         ++z)
    {
        const double coarseZ =
            static_cast<double>(z)
            *
            static_cast<double>(coarseNz - 1)
            /
            static_cast<double>(fineNz - 1);

        for (int y = 1;
             y < fineNy - 1;
             ++y)
        {
            const double coarseY =
                static_cast<double>(y)
                *
                static_cast<double>(coarseNy - 1)
                /
                static_cast<double>(fineNy - 1);

            for (int x = 1;
                 x < fineNx - 1;
                 ++x)
            {
                const double coarseX =
                    static_cast<double>(x)
                    *
                    static_cast<double>(coarseNx - 1)
                    /
                    static_cast<double>(fineNx - 1);

                const std::size_t fineIndex =
                    static_cast<std::size_t>(x)
                    +
                    static_cast<std::size_t>(fineNx)
                    *
                    (
                        static_cast<std::size_t>(y)
                        +
                        static_cast<std::size_t>(fineNy)
                        *
                        static_cast<std::size_t>(z)
                    );

                fineSolution[fineIndex] +=
                    sampleTrilinear(
                        coarseCorrection,
                        coarseNx,
                        coarseNy,
                        coarseNz,
                        coarseX,
                        coarseY,
                        coarseZ
                    );
            }
        }
    }
}


// ================================================================
// INTERPOLACIÓN TRILINEAL
// ================================================================

double DFTPotential::sampleTrilinear(
    const std::vector<double>& field,
    int nx,
    int ny,
    int nz,
    double x,
    double y,
    double z
) const
{
    x =
        std::clamp(
            x,
            0.0,
            static_cast<double>(nx - 1)
        );

    y =
        std::clamp(
            y,
            0.0,
            static_cast<double>(ny - 1)
        );

    z =
        std::clamp(
            z,
            0.0,
            static_cast<double>(nz - 1)
        );

    const int x0 =
        static_cast<int>(
            std::floor(x)
        );

    const int y0 =
        static_cast<int>(
            std::floor(y)
        );

    const int z0 =
        static_cast<int>(
            std::floor(z)
        );

    const int x1 =
        std::min(
            x0 + 1,
            nx - 1
        );

    const int y1 =
        std::min(
            y0 + 1,
            ny - 1
        );

    const int z1 =
        std::min(
            z0 + 1,
            nz - 1
        );

    const double tx =
        x - static_cast<double>(x0);

    const double ty =
        y - static_cast<double>(y0);

    const double tz =
        z - static_cast<double>(z0);

    auto index =
        [nx, ny](int ix, int iy, int iz)
        {
            return
                static_cast<std::size_t>(ix)
                +
                static_cast<std::size_t>(nx)
                *
                (
                    static_cast<std::size_t>(iy)
                    +
                    static_cast<std::size_t>(ny)
                    *
                    static_cast<std::size_t>(iz)
                );
        };

    const double c000 =
        field[index(x0, y0, z0)];

    const double c100 =
        field[index(x1, y0, z0)];

    const double c010 =
        field[index(x0, y1, z0)];

    const double c110 =
        field[index(x1, y1, z0)];

    const double c001 =
        field[index(x0, y0, z1)];

    const double c101 =
        field[index(x1, y0, z1)];

    const double c011 =
        field[index(x0, y1, z1)];

    const double c111 =
        field[index(x1, y1, z1)];

    const double c00 =
        c000 * (1.0 - tx) +
        c100 * tx;

    const double c10 =
        c010 * (1.0 - tx) +
        c110 * tx;

    const double c01 =
        c001 * (1.0 - tx) +
        c101 * tx;

    const double c11 =
        c011 * (1.0 - tx) +
        c111 * tx;

    const double c0 =
        c00 * (1.0 - ty) +
        c10 * ty;

    const double c1 =
        c01 * (1.0 - ty) +
        c11 * ty;

    return
        c0 * (1.0 - tz) +
        c1 * tz;
}


// ================================================================
// V-CYCLE MULTIGRID
// ================================================================

void DFTPotential::vCycle(
    std::vector<double>& solution,
    const std::vector<double>& rhs,
    int nx,
    int ny,
    int nz,
    double hx,
    double hy,
    double hz
) const
{
    /*
        En el nivel más grueso resolvemos directamente
        mediante suficientes barridos de Gauss-Seidel.

        La solución del problema de corrección tiene
        condiciones de frontera homogéneas:

            e = 0

        por lo que los valores de frontera de
        "solution" permanecen en cero.
    */

    if (nx <= 3 ||
        ny <= 3 ||
        nz <= 3)
    {
        smooth(
            solution,
            rhs,
            nx,
            ny,
            nz,
            hx,
            hy,
            hz,
            COARSE_SMOOTH
        );

        return;
    }

    // ------------------------------------------------------------
    // 1. Pre-smoothing
    // ------------------------------------------------------------

    smooth(
        solution,
        rhs,
        nx,
        ny,
        nz,
        hx,
        hy,
        hz,
        PRE_SMOOTH
    );

    // ------------------------------------------------------------
    // 2. Residuo
    // ------------------------------------------------------------

    const std::vector<double> residual =
        calculateResidual(
            solution,
            rhs,
            nx,
            ny,
            nz,
            hx,
            hy,
            hz
        );

    // ------------------------------------------------------------
    // 3. Construcción del grid grueso
    //
    // Se conservan los extremos físicos del dominio.
    // ------------------------------------------------------------

    const int coarseNx =
        nx / 2 + 1;

    const int coarseNy =
        ny / 2 + 1;

    const int coarseNz =
        nz / 2 + 1;

    const double domainX =
        static_cast<double>(nx - 1) *
        hx;

    const double domainY =
        static_cast<double>(ny - 1) *
        hy;

    const double domainZ =
        static_cast<double>(nz - 1) *
        hz;

    const double coarseHx =
        domainX /
        static_cast<double>(coarseNx - 1);

    const double coarseHy =
        domainY /
        static_cast<double>(coarseNy - 1);

    const double coarseHz =
        domainZ /
        static_cast<double>(coarseNz - 1);

    // ------------------------------------------------------------
    // 4. Restriction del residuo
    // ------------------------------------------------------------

    const std::vector<double> coarseRhs =
        restrictResidual(
            residual,
            nx,
            ny,
            nz,
            coarseNx,
            coarseNy,
            coarseNz
        );

    // ------------------------------------------------------------
    // 5. Problema de corrección en el grid grueso
    // ------------------------------------------------------------

    std::vector<double> coarseCorrection(
        static_cast<std::size_t>(coarseNx)
        *
        static_cast<std::size_t>(coarseNy)
        *
        static_cast<std::size_t>(coarseNz),
        0.0
    );

    vCycle(
        coarseCorrection,
        coarseRhs,
        coarseNx,
        coarseNy,
        coarseNz,
        coarseHx,
        coarseHy,
        coarseHz
    );

    // ------------------------------------------------------------
    // 6. Prolongation + corrección
    // ------------------------------------------------------------

    prolongateAndAdd(
        coarseCorrection,
        coarseNx,
        coarseNy,
        coarseNz,
        solution,
        nx,
        ny,
        nz
    );

    // ------------------------------------------------------------
    // 7. Post-smoothing
    // ------------------------------------------------------------

    smooth(
        solution,
        rhs,
        nx,
        ny,
        nz,
        hx,
        hy,
        hz,
        POST_SMOOTH
    );
}


// ================================================================
// SOLVER POISSON
// ================================================================

void DFTPotential::solvePoisson(
    const DFTDensity& density,
    std::vector<double>& potential,
    int maxCycles
)
{
    const int nx =
        grid->getNx();

    const int ny =
        grid->getNy();

    const int nz =
        grid->getNz();

    const double h =
        grid->getSpacing();

    const std::size_t count =
        grid->getPointCount();

    if (potential.size() != count)
    {
        potential.assign(
            count,
            0.0
        );
    }

    /*
        El potencial anterior se reutiliza como
        solución inicial.

        Esto es especialmente importante dentro
        del SCF, donde rho[n] normalmente es
        cercana a rho[n-1].
    */

    applyBoundaryConditions(
        potential,
        density
    );

    /*
        Construimos:

            ∇²V = -4πρ

        por tanto:

            rhs = -4πρ
    */

    std::vector<double> rhs(
        count,
        0.0
    );

    for (std::size_t i = 0;
         i < count;
         ++i)
    {
        const double rho =
            density.get(i);

        if (!std::isfinite(rho))
        {
            throw std::runtime_error(
                "DFTPotential: density contains a non-finite value."
            );
        }

        rhs[i] =
            -4.0 *
            PI *
            rho;
    }

    /*
        Las condiciones de frontera físicas no forman
        parte del residuo interior.
    */

    const std::vector<double> initialResidual =
        calculateResidual(
            potential,
            rhs,
            nx,
            ny,
            nz,
            h,
            h,
            h
        );

    const double initialNorm =
        calculateResidualNorm(
            initialResidual,
            nx,
            ny,
            nz,
            h,
            h,
            h
        );

    poissonResidual =
        initialNorm;

    poissonCycles =
        0;

    poissonConverged =
        false;

    /*
        Si ya estamos dentro de la tolerancia,
        no necesitamos ejecutar ningún V-cycle.
    */

    if (initialNorm <= POISSON_TOLERANCE)
    {
        poissonConverged =
            true;

        return;
    }

    /*
        Utilizamos simultáneamente:

            residuo absoluto

        y

            residuo relativo.

        La condición:

            residual <= tolerance * max(1, initialNorm)

        evita problemas cuando el residuo inicial
        es extremadamente pequeño.
    */

    const double convergenceThreshold =
        POISSON_TOLERANCE *
        std::max(
            1.0,
            initialNorm
        );

    for (int cycle = 1;
         cycle <= maxCycles;
         ++cycle)
    {
        vCycle(
            potential,
            rhs,
            nx,
            ny,
            nz,
            h,
            h,
            h
        );

        /*
            Las condiciones físicas de frontera
            pertenecen únicamente al problema
            original de Poisson.

            El V-cycle resuelve las correcciones
            con frontera homogénea, por lo que
            volvemos a imponer la frontera física.
        */

        applyBoundaryConditions(
            potential,
            density
        );

        const std::vector<double> residual =
            calculateResidual(
                potential,
                rhs,
                nx,
                ny,
                nz,
                h,
                h,
                h
            );

        poissonResidual =
            calculateResidualNorm(
                residual,
                nx,
                ny,
                nz,
                h,
                h,
                h
            );

        poissonCycles =
            cycle;

        if (poissonResidual <=
            convergenceThreshold)
        {
            poissonConverged =
                true;

            return;
        }
    }

    /*
        Llegamos al máximo de V-cycles sin
        alcanzar la tolerancia.

        No lanzamos una excepción aquí:
        el potencial obtenido sigue siendo la
        mejor aproximación disponible y el estado
        de convergencia queda explícitamente
        registrado.
    */
}


// ================================================================
// INTERCAMBIO LDA
// ================================================================

void DFTPotential::calculateExchangePotential(
    const DFTDensity& density
)
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "DFTPotential: not initialized."
        );
    }

    if (density.getGrid().getPointCount() !=
        grid->getPointCount())
    {
        throw std::invalid_argument(
            "DFTPotential: incompatible density grid."
        );
    }

    const double coefficient =
        -std::cbrt(
            3.0 / PI
        );

    const std::size_t count =
        grid->getPointCount();

    for (std::size_t i = 0;
         i < count;
         ++i)
    {
        const double rho =
            std::max(
                density.get(i),
                0.0
            );

        if (rho <= MIN_DENSITY)
        {
            exchangePotential[i] =
                0.0;
        }
        else
        {
            exchangePotential[i] =
                coefficient *
                std::cbrt(rho);
        }
    }
}


// ================================================================
// ACCESO
// ================================================================

double DFTPotential::getExternal(
    std::size_t index
) const
{
    return externalPotential.at(
        index
    );
}


double DFTPotential::getHartree(
    std::size_t index
) const
{
    return hartreePotential.at(
        index
    );
}


double DFTPotential::getExchange(
    std::size_t index
) const
{
    return exchangePotential.at(
        index
    );
}


double DFTPotential::getTotalElectronic(
    std::size_t index
) const
{
    return
        externalPotential.at(index)
        +
        hartreePotential.at(index)
        +
        exchangePotential.at(index);
}


const std::vector<double>&
DFTPotential::getExternalPotential() const
{
    return externalPotential;
}


const std::vector<double>&
DFTPotential::getHartreePotential() const
{
    return hartreePotential;
}


const std::vector<double>&
DFTPotential::getExchangePotential() const
{
    return exchangePotential;
}


// ================================================================
// DIAGNÓSTICO POISSON
// ================================================================

double DFTPotential::getPoissonResidual() const
{
    return poissonResidual;
}


int DFTPotential::getPoissonCycles() const
{
    return poissonCycles;
}


bool DFTPotential::hasPoissonConverged() const
{
    return poissonConverged;
}


// ================================================================
// ENERGÍA ELECTRÓN - NÚCLEO
// ================================================================

double DFTPotential::calculateElectronNuclearEnergy(
    const DFTDensity& density
) const
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "DFTPotential: not initialized."
        );
    }

    const double dV =
        grid->getVolumeElement();

    double energy =
        0.0;

    for (std::size_t i = 0;
         i < grid->getPointCount();
         ++i)
    {
        energy +=
            density.get(i) *
            externalPotential[i] *
            dV;
    }

    return energy;
}


// ================================================================
// ENERGÍA HARTREE
// ================================================================

double DFTPotential::calculateHartreeEnergy(
    const DFTDensity& density
) const
{
    const double dV =
        grid->getVolumeElement();

    double energy =
        0.0;

    for (std::size_t i = 0;
         i < grid->getPointCount();
         ++i)
    {
        energy +=
            density.get(i) *
            hartreePotential[i] *
            dV;
    }

    return 0.5 * energy;
}


// ================================================================
// ENERGÍA DE INTERCAMBIO
// ================================================================

double DFTPotential::calculateExchangeEnergy(
    const DFTDensity& density
) const
{
    const double coefficient =
        -0.75 *
        std::cbrt(
            3.0 / PI
        );

    const double dV =
        grid->getVolumeElement();

    double energy =
        0.0;

    for (std::size_t i = 0;
         i < grid->getPointCount();
         ++i)
    {
        const double rho =
            std::max(
                density.get(i),
                0.0
            );

        if (rho > MIN_DENSITY)
        {
            energy +=
                coefficient *
                std::pow(
                    rho,
                    4.0 / 3.0
                ) *
                dV;
        }
    }

    return energy;
}


// ================================================================
// REPULSIÓN NÚCLEO - NÚCLEO
// ================================================================

double DFTPotential::calculateNuclearRepulsionEnergy(
    const std::vector<int>& nuclearCharges,
    const std::vector<glm::dvec3>& nuclearPositions
) const
{
    if (nuclearCharges.size() !=
        nuclearPositions.size())
    {
        throw std::invalid_argument(
            "DFTPotential: nuclear charges and positions have different sizes."
        );
    }

    double energy =
        0.0;

    for (std::size_t i = 0;
         i < nuclearCharges.size();
         ++i)
    {
        if (nuclearCharges[i] <= 0)
        {
            throw std::invalid_argument(
                "DFTPotential: nuclear charge must be positive."
            );
        }

        for (std::size_t j = i + 1;
             j < nuclearCharges.size();
             ++j)
        {
            const glm::dvec3 difference =
                nuclearPositions[i] -
                nuclearPositions[j];

            const double distance =
                glm::length(
                    difference
                );

            if (distance <= MIN_DISTANCE)
            {
                throw std::runtime_error(
                    "DFTPotential: nuclei overlap."
                );
            }

            energy +=
                static_cast<double>(
                    nuclearCharges[i]
                )
                *
                static_cast<double>(
                    nuclearCharges[j]
                )
                /
                distance;
        }
    }

    return energy;
}


// ================================================================
// FUERZAS
// ================================================================

std::vector<glm::dvec3>
DFTPotential::calculateNuclearForces(
    const DFTDensity& density,
    const std::vector<int>& nuclearCharges,
    const std::vector<glm::dvec3>& nuclearPositions
) const
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "DFTPotential: not initialized."
        );
    }

    if (nuclearCharges.size() !=
        nuclearPositions.size())
    {
        throw std::invalid_argument(
            "DFTPotential: nuclear charges and positions have different sizes."
        );
    }

    std::vector<glm::dvec3> forces(
        nuclearPositions.size(),
        glm::dvec3(0.0)
    );

    const double dV =
        grid->getVolumeElement();

    /*
        Fuerza electrón-núcleo consistente
        con el potencial suavizado actual.
    */

    for (std::size_t nucleus = 0;
         nucleus < nuclearPositions.size();
         ++nucleus)
    {
        const double Z =
            static_cast<double>(
                nuclearCharges[nucleus]
            );

        glm::dvec3 electronicForce(
            0.0
        );

        for (std::size_t i = 0;
             i < grid->getPointCount();
             ++i)
        {
            const glm::dvec3 difference =
                grid->getPosition(i) -
                nuclearPositions[nucleus];

            const double r2 =
                glm::dot(
                    difference,
                    difference
                );

            const double denominator =
                std::pow(
                    r2 +
                    softening * softening,
                    1.5
                );

            electronicForce +=
                Z *
                density.get(i) *
                difference /
                denominator *
                dV;
        }

        forces[nucleus] -=
            electronicForce;
    }

    /*
        Fuerza núcleo-núcleo.
    */

    for (std::size_t i = 0;
         i < nuclearPositions.size();
         ++i)
    {
        for (std::size_t j = i + 1;
             j < nuclearPositions.size();
             ++j)
        {
            const glm::dvec3 difference =
                nuclearPositions[i] -
                nuclearPositions[j];

            const double r2 =
                glm::dot(
                    difference,
                    difference
                );

            const double distance =
                std::sqrt(r2);

            if (distance <= MIN_DISTANCE)
            {
                throw std::runtime_error(
                    "DFTPotential: nuclei overlap."
                );
            }

            const glm::dvec3 force =
                static_cast<double>(
                    nuclearCharges[i]
                )
                *
                static_cast<double>(
                    nuclearCharges[j]
                )
                *
                difference /
                (
                    distance *
                    distance *
                    distance
                );

            forces[i] +=
                force;

            forces[j] -=
                force;
        }
    }

    return forces;
}


// ================================================================
// SOFTENING
// ================================================================

void DFTPotential::setSoftening(
    double value
)
{
    if (!std::isfinite(value) ||
        value <= 0.0)
    {
        throw std::invalid_argument(
            "DFTPotential: softening must be positive and finite."
        );
    }

    softening =
        value;
}


double DFTPotential::getSoftening() const
{
    return softening;
}


// ================================================================
// GRID
// ================================================================

const DFTGrid&
DFTPotential::getGrid() const
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "DFTPotential: not initialized."
        );
    }

    return *grid;
}