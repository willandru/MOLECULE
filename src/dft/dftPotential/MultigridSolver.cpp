#include "MultigridSolver.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace
{
    constexpr double PI =
        3.1415926535897932384626433832795;

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

MultigridSolver::MultigridSolver()
    :
    grid(nullptr),
    residual(
        std::numeric_limits<double>::infinity()
    ),
    cycles(0),
    converged(false)
{
}


MultigridSolver::MultigridSolver(
    const DFTGrid& grid
)
    :
    grid(nullptr),
    residual(
        std::numeric_limits<double>::infinity()
    ),
    cycles(0),
    converged(false)
{
    initialize(grid);
}


// ================================================================
// INICIALIZACIÓN
// ================================================================

void MultigridSolver::initialize(
    const DFTGrid& grid
)
{
    if (grid.getPointCount() == 0)
    {
        throw std::invalid_argument(
            "MultigridSolver: grid has no points."
        );
    }

    if (grid.getNx() < 3 ||
        grid.getNy() < 3 ||
        grid.getNz() < 3)
    {
        throw std::invalid_argument(
            "MultigridSolver: grid dimensions must be at least 3 in every direction."
        );
    }

    if (!std::isfinite(grid.getSpacing()) ||
        grid.getSpacing() <= 0.0)
    {
        throw std::invalid_argument(
            "MultigridSolver: grid spacing must be positive and finite."
        );
    }

    this->grid =
        &grid;

    residual =
        std::numeric_limits<double>::infinity();

    cycles =
        0;

    converged =
        false;
}


// ================================================================
// SUAVIZADO GAUSS-SEIDEL
// ================================================================

void MultigridSolver::smooth(
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
        2.0 *
        (
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
MultigridSolver::calculateResidual(
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
            "MultigridSolver: invalid Poisson vector size."
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

double MultigridSolver::calculateResidualNorm(
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
            "MultigridSolver: invalid residual size."
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
// INTERPOLACIÓN TRILINEAL
// ================================================================

double MultigridSolver::sampleTrilinear(
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
        x -
        static_cast<double>(x0);

    const double ty =
        y -
        static_cast<double>(y0);

    const double tz =
        z -
        static_cast<double>(z0);

    auto index =
        [nx, ny](
            int ix,
            int iy,
            int iz
        )
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
// RESTRICCIÓN
// ================================================================

std::vector<double>
MultigridSolver::restrictResidual(
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
            "MultigridSolver: invalid fine residual size."
        );
    }

    const std::size_t coarseCount =
        static_cast<std::size_t>(coarseNx)
        *
        static_cast<std::size_t>(coarseNy)
        *
        static_cast<std::size_t>(coarseNz);

    std::vector<double> coarseResidual(
        coarseCount,
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
// PROLONGACIÓN
// ================================================================

void MultigridSolver::prolongateAndAdd(
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
            "MultigridSolver: invalid multigrid vector size."
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
// V-CYCLE
// ================================================================

void MultigridSolver::vCycle(
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
        El problema de corrección tiene condiciones
        de frontera homogéneas.

        Por eso este método nunca modifica los
        valores de frontera.
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
    // 3. Grid grueso
    //
    // Se conservan los extremos físicos.
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
    // 4. Restricción
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
    // 5. Problema de corrección
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
    // 6. Prolongación y corrección
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
// SOLVER
// ================================================================

void MultigridSolver::solve(
    const DFTDensity& density,
    std::vector<double>& potential,
    int maxCycles
)
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "MultigridSolver: not initialized."
        );
    }

    if (maxCycles <= 0)
    {
        throw std::invalid_argument(
            "MultigridSolver: maximum number of cycles must be positive."
        );
    }

    if (density.getGrid().getPointCount() !=
        grid->getPointCount())
    {
        throw std::invalid_argument(
            "MultigridSolver: incompatible density grid."
        );
    }

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

    // ------------------------------------------------------------
    // Construcción del lado derecho
    //
    // ∇²V_H = -4πρ
    //
    // A V = rhs
    // ------------------------------------------------------------

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
                "MultigridSolver: density contains a non-finite value."
            );
        }

        rhs[i] =
            -4.0 *
            PI *
            rho;
    }

    // ------------------------------------------------------------
    // Residuo inicial
    // ------------------------------------------------------------

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

    residual =
        initialNorm;

    cycles =
        0;

    converged =
        false;

    if (initialNorm <= POISSON_TOLERANCE)
    {
        converged =
            true;

        return;
    }

    const double convergenceThreshold =
        POISSON_TOLERANCE *
        std::max(
            1.0,
            initialNorm
        );

    // ------------------------------------------------------------
    // V-cycles
    // ------------------------------------------------------------

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

        const std::vector<double> currentResidual =
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

        residual =
            calculateResidualNorm(
                currentResidual,
                nx,
                ny,
                nz,
                h,
                h,
                h
            );

        cycles =
            cycle;

        if (residual <=
            convergenceThreshold)
        {
            converged =
                true;

            return;
        }
    }
}


// ================================================================
// DIAGNÓSTICO
// ================================================================

double MultigridSolver::getResidual() const
{
    return residual;
}


int MultigridSolver::getCycles() const
{
    return cycles;
}


bool MultigridSolver::hasConverged() const
{
    return converged;
}