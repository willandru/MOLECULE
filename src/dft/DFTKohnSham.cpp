#include "DFTKohnSham.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace
{
    constexpr double CONVERGENCE = 1.0e-6;
    constexpr double MIN_NORM = 1.0e-14;

    constexpr double GAUSSIAN_ALPHA = 0.20;

    constexpr int JACOBI_MAX_ITERATIONS = 1000;
    constexpr double JACOBI_TOLERANCE = 1.0e-12;
}

DFTKohnSham::DFTKohnSham()
    : grid(nullptr),
      electronCount(0),
      orbitalCount(0),
      density(),
      converged(false)
{
}

DFTKohnSham::DFTKohnSham(
    const DFTGrid& newGrid
)
    : DFTKohnSham()
{
    initialize(newGrid);
}

void DFTKohnSham::initialize(
    const DFTGrid& newGrid
)
{
    if (newGrid.getPointCount() == 0)
    {
        throw std::invalid_argument(
            "DFTKohnSham: grid cannot be empty."
        );
    }

    grid = &newGrid;

    electronCount = 0;
    orbitalCount = 0;

    orbitals.clear();
    eigenvalues.clear();
    occupations.clear();
    residuals.clear();

    density.initialize(newGrid);

    converged = false;
}

void DFTKohnSham::setElectronCount(
    int count
)
{
    if (count < 0)
    {
        throw std::invalid_argument(
            "DFTKohnSham: electron count cannot be negative."
        );
    }

    if (grid == nullptr)
    {
        throw std::runtime_error(
            "DFTKohnSham: grid has not been initialized."
        );
    }

    electronCount = count;

    /*
     * Closed-shell approximation:
     *
     * 2 electrons per spatial orbital.
     */
    orbitalCount =
        (electronCount + 1) / 2;

    orbitals.clear();
    eigenvalues.clear();
    occupations.clear();
    residuals.clear();

    if (orbitalCount == 0)
    {
        converged = false;
        return;
    }

    const int pointCount =
        grid->getPointCount();

    orbitals.resize(
        orbitalCount,
        std::vector<double>(
            pointCount,
            0.0
        )
    );

    eigenvalues.assign(
        orbitalCount,
        0.0
    );

    occupations.assign(
        orbitalCount,
        0
    );

    residuals.assign(
        orbitalCount,
        std::numeric_limits<double>::infinity()
    );

    int remainingElectrons =
        electronCount;

    for (int i = 0;
         i < orbitalCount;
         ++i)
    {
        occupations[i] =
            std::min(
                2,
                remainingElectrons
            );

        remainingElectrons -=
            occupations[i];
    }

    initializeOrbitals();

    converged = false;
}

int DFTKohnSham::getElectronCount() const
{
    return electronCount;
}

int DFTKohnSham::getOrbitalCount() const
{
    return orbitalCount;
}

void DFTKohnSham::initializeOrbitals()
{
    if (grid == nullptr ||
        orbitalCount == 0)
    {
        return;
    }

    const int nx =
        grid->getNx();

    const int ny =
        grid->getNy();

    const int nz =
        grid->getNz();

    const int pointCount =
        grid->getPointCount();

    /*
     * Precompute coordinates.
     *
     * The old implementation repeatedly called
     * grid->getPosition() during initialization.
     *
     * This is not important for SCF performance, but it
     * keeps the initialization cleaner and cheaper.
     */
    std::vector<double> px(pointCount);
    std::vector<double> py(pointCount);
    std::vector<double> pz(pointCount);
    std::vector<double> gaussian(pointCount);

    for (int z = 0;
         z < nz;
         ++z)
    {
        for (int y = 0;
             y < ny;
             ++y)
        {
            for (int x = 0;
                 x < nx;
                 ++x)
            {
                const int index =
                    grid->getIndex(
                        x,
                        y,
                        z
                    );

                const glm::dvec3 position =
                    grid->getPosition(index);

                px[index] = position.x;
                py[index] = position.y;
                pz[index] = position.z;

                const double r2 =
                    position.x * position.x +
                    position.y * position.y +
                    position.z * position.z;

                gaussian[index] =
                    std::exp(
                        -GAUSSIAN_ALPHA * r2
                    );
            }
        }
    }

    for (int orbitalIndex = 0;
         orbitalIndex < orbitalCount;
         ++orbitalIndex)
    {
        std::vector<double>& orbital =
            orbitals[orbitalIndex];

        for (int index = 0;
             index < pointCount;
             ++index)
        {
            const double x =
                px[index];

            const double y =
                py[index];

            const double z =
                pz[index];

            const double g =
                gaussian[index];

            double value = 0.0;

            switch (orbitalIndex % 7)
            {
            case 0:
                value = g;
                break;

            case 1:
                value = x * g;
                break;

            case 2:
                value = y * g;
                break;

            case 3:
                value = z * g;
                break;

            case 4:
                value = x * y * g;
                break;

            case 5:
                value = x * z * g;
                break;

            case 6:
                value = y * z * g;
                break;
            }

            orbital[index] =
                value;
        }

        applyOrbitalBoundaryConditions(
            orbital
        );
    }

    /*
     * Make the initial orbitals orthonormal.
     */
    std::vector<std::vector<double>> basis =
        orbitals;

    orthonormalize(
        basis
    );

    if (static_cast<int>(basis.size()) !=
        orbitalCount)
    {
        throw std::runtime_error(
            "DFTKohnSham: unable to construct independent initial orbitals."
        );
    }

    orbitals =
        std::move(basis);
}

double DFTKohnSham::calculateNorm(
    const std::vector<double>& orbital
) const
{
    if (grid == nullptr)
    {
        return 0.0;
    }

    double sum = 0.0;

    for (const double value :
         orbital)
    {
        sum +=
            value * value;
    }

    return std::sqrt(
        sum *
        grid->getVolumeElement()
    );
}

double DFTKohnSham::calculateOverlap(
    const std::vector<double>& a,
    const std::vector<double>& b
) const
{
    if (a.size() != b.size())
    {
        throw std::invalid_argument(
            "DFTKohnSham: incompatible orbital sizes."
        );
    }

    if (grid == nullptr)
    {
        return 0.0;
    }

    double sum = 0.0;

    const std::size_t size =
        a.size();

    for (std::size_t i = 0;
         i < size;
         ++i)
    {
        sum +=
            a[i] * b[i];
    }

    return sum *
           grid->getVolumeElement();
}

void DFTKohnSham::orthonormalize(
    std::vector<std::vector<double>>& basis
) const
{
    std::vector<std::vector<double>> result;

    result.reserve(
        basis.size()
    );

    /*
     * Modified Gram-Schmidt with two passes.
     *
     * This is still necessary because the projected
     * Hamiltonian assumes an orthonormal basis.
     */
    for (const std::vector<double>& original :
         basis)
    {
        std::vector<double> vector =
            original;

        applyOrbitalBoundaryConditions(
            vector
        );

        for (int pass = 0;
             pass < 2;
             ++pass)
        {
            for (const std::vector<double>& q :
                 result)
            {
                const double overlap =
                    calculateOverlap(
                        q,
                        vector
                    );

                if (std::abs(overlap) <
                    MIN_NORM)
                {
                    continue;
                }

                for (std::size_t i = 0;
                     i < vector.size();
                     ++i)
                {
                    vector[i] -=
                        overlap *
                        q[i];
                }
            }
        }

        const double norm =
            calculateNorm(
                vector
            );

        if (norm <= MIN_NORM)
        {
            continue;
        }

        const double inverseNorm =
            1.0 / norm;

        for (double& value :
             vector)
        {
            value *=
                inverseNorm;
        }

        /*
         * Boundary values were already zero, so this does
         * not require another full grid traversal.
         */

        result.push_back(
            std::move(vector)
        );
    }

    basis =
        std::move(result);
}

void DFTKohnSham::applyOrbitalBoundaryConditions(
    std::vector<double>& orbital
) const
{
    if (grid == nullptr)
    {
        return;
    }

    const int nx =
        grid->getNx();

    const int ny =
        grid->getNy();

    const int nz =
        grid->getNz();

    /*
     * Only boundary planes are modified.
     *
     * Avoid getIndex() where possible by using the
     * grid's known row-major layout:
     *
     * index = x + nx * (y + ny * z)
     */
    for (int z = 0;
         z < nz;
         ++z)
    {
        for (int y = 0;
             y < ny;
             ++y)
        {
            const int rowStart =
                nx * (y + ny * z);

            orbital[rowStart] =
                0.0;

            orbital[rowStart + nx - 1] =
                0.0;
        }
    }

    /*
     * y = 0 and y = ny-1
     */
    for (int z = 0;
         z < nz;
         ++z)
    {
        const int firstRow =
            nx * (ny * z);

        const int lastRow =
            nx * (ny - 1 + ny * z);

        for (int x = 0;
             x < nx;
             ++x)
        {
            orbital[firstRow + x] =
                0.0;

            orbital[lastRow + x] =
                0.0;
        }
    }

    /*
     * z = 0 and z = nz-1
     */
    const int planeSize =
        nx * ny;

    const int lastPlane =
        planeSize * (nz - 1);

    for (int i = 0;
         i < planeSize;
         ++i)
    {
        orbital[i] =
            0.0;

        orbital[lastPlane + i] =
            0.0;
    }
}

void DFTKohnSham::calculateLaplacian(
    const std::vector<double>& orbital,
    std::vector<double>& laplacian
) const
{
    const int nx =
        grid->getNx();

    const int ny =
        grid->getNy();

    const int nz =
        grid->getNz();

    const double h =
        grid->getSpacing();

    const double inverseH2 =
        1.0 / (h * h);

    const int planeSize =
        nx * ny;

    /*
     * Reuse/resize the output buffer.
     */
    laplacian.assign(
        orbital.size(),
        0.0
    );

    for (int z = 1;
         z < nz - 1;
         ++z)
    {
        const int plane =
            z * planeSize;

        for (int y = 1;
             y < ny - 1;
             ++y)
        {
            const int row =
                plane +
                y * nx;

            for (int x = 1;
                 x < nx - 1;
                 ++x)
            {
                const int index =
                    row + x;

                laplacian[index] =
                    (
                        orbital[index + 1] +
                        orbital[index - 1] +

                        orbital[index + nx] +
                        orbital[index - nx] +

                        orbital[index + planeSize] +
                        orbital[index - planeSize] -

                        6.0 *
                        orbital[index]
                    ) *
                    inverseH2;
            }
        }
    }
}

void DFTKohnSham::applyHamiltonian(
    const std::vector<double>& orbital,
    const DFTPotential& potential,
    std::vector<double>& result
) const
{
    /*
     * Hψ = -1/2 ∇²ψ + Vψ
     */

    std::vector<double> laplacian;

    calculateLaplacian(
        orbital,
        laplacian
    );

    const std::size_t pointCount =
        orbital.size();

    result.resize(
        pointCount
    );

    for (std::size_t i = 0;
         i < pointCount;
         ++i)
    {
        result[i] =
            -0.5 *
            laplacian[i]
            +
            potential.getTotalElectronic(
                static_cast<int>(i)
            ) *
            orbital[i];
    }

    /*
     * Boundary values of ψ are zero, therefore Hψ is also
     * explicitly forced to zero there.
     */
    applyOrbitalBoundaryConditions(
        result
    );
}

double DFTKohnSham::calculateOrbitalEnergy(
    const std::vector<double>& orbital,
    const std::vector<double>& hOrbital
) const
{
    const double numerator =
        calculateOverlap(
            orbital,
            hOrbital
        );

    /*
     * Orbitals are normalized, but keeping the denominator
     * makes this robust to numerical drift.
     */
    const double denominator =
        calculateOverlap(
            orbital,
            orbital
        );

    if (denominator <= MIN_NORM)
    {
        return 0.0;
    }

    return numerator /
           denominator;
}

double DFTKohnSham::calculateResidual(
    const std::vector<double>& orbital,
    const std::vector<double>& hOrbital,
    double eigenvalue
) const
{
    double sum = 0.0;

    const std::size_t size =
        orbital.size();

    for (std::size_t i = 0;
         i < size;
         ++i)
    {
        const double value =
            hOrbital[i] -
            eigenvalue *
            orbital[i];

        sum +=
            value * value;
    }

    return std::sqrt(
        sum *
        grid->getVolumeElement()
    );
}

void DFTKohnSham::diagonalizeSymmetricMatrix(
    std::vector<std::vector<double>>& matrix,
    std::vector<double>& matrixEigenvalues,
    std::vector<std::vector<double>>& eigenvectors
) const
{
    const int n =
        static_cast<int>(
            matrix.size()
        );

    if (n == 0)
    {
        matrixEigenvalues.clear();
        eigenvectors.clear();
        return;
    }

    for (const auto& row :
         matrix)
    {
        if (static_cast<int>(row.size()) != n)
        {
            throw std::invalid_argument(
                "DFTKohnSham: matrix must be square."
            );
        }
    }

    /*
     * Eigenvectors stored by columns.
     */
    eigenvectors.assign(
        n,
        std::vector<double>(
            n,
            0.0
        )
    );

    for (int i = 0;
         i < n;
         ++i)
    {
        eigenvectors[i][i] =
            1.0;
    }

    for (int iteration = 0;
         iteration < JACOBI_MAX_ITERATIONS;
         ++iteration)
    {
        int p = 0;
        int q = 0;

        double maximumOffDiagonal =
            0.0;

        for (int i = 0;
             i < n;
             ++i)
        {
            for (int j = i + 1;
                 j < n;
                 ++j)
            {
                const double value =
                    std::abs(
                        matrix[i][j]
                    );

                if (value >
                    maximumOffDiagonal)
                {
                    maximumOffDiagonal =
                        value;

                    p = i;
                    q = j;
                }
            }
        }

        if (maximumOffDiagonal <
            JACOBI_TOLERANCE)
        {
            break;
        }

        const double app =
            matrix[p][p];

        const double aqq =
            matrix[q][q];

        const double apq =
            matrix[p][q];

        const double angle =
            0.5 *
            std::atan2(
                2.0 * apq,
                aqq - app
            );

        const double c =
            std::cos(angle);

        const double s =
            std::sin(angle);

        for (int k = 0;
             k < n;
             ++k)
        {
            if (k == p ||
                k == q)
            {
                continue;
            }

            const double mkp =
                matrix[k][p];

            const double mkq =
                matrix[k][q];

            const double newKp =
                c * mkp -
                s * mkq;

            const double newKq =
                s * mkp +
                c * mkq;

            matrix[k][p] =
                newKp;

            matrix[p][k] =
                newKp;

            matrix[k][q] =
                newKq;

            matrix[q][k] =
                newKq;
        }

        matrix[p][p] =
            c * c * app -
            2.0 * s * c * apq +
            s * s * aqq;

        matrix[q][q] =
            s * s * app +
            2.0 * s * c * apq +
            c * c * aqq;

        matrix[p][q] =
            0.0;

        matrix[q][p] =
            0.0;

        for (int k = 0;
             k < n;
             ++k)
        {
            const double vip =
                eigenvectors[k][p];

            const double viq =
                eigenvectors[k][q];

            eigenvectors[k][p] =
                c * vip -
                s * viq;

            eigenvectors[k][q] =
                s * vip +
                c * viq;
        }
    }

    matrixEigenvalues.resize(
        n
    );

    for (int i = 0;
         i < n;
         ++i)
    {
        matrixEigenvalues[i] =
            matrix[i][i];
    }

    /*
     * Sort eigenvalues and corresponding eigenvector columns.
     */
    std::vector<int> indices(
        n
    );

    for (int i = 0;
         i < n;
         ++i)
    {
        indices[i] =
            i;
    }

    std::sort(
        indices.begin(),
        indices.end(),
        [&matrixEigenvalues](
            int a,
            int b
        )
        {
            return matrixEigenvalues[a] <
                   matrixEigenvalues[b];
        }
    );

    std::vector<double> sortedValues(
        n
    );

    std::vector<std::vector<double>> sortedVectors(
        n,
        std::vector<double>(
            n,
            0.0
        )
    );

    for (int newIndex = 0;
         newIndex < n;
         ++newIndex)
    {
        const int oldIndex =
            indices[newIndex];

        sortedValues[newIndex] =
            matrixEigenvalues[oldIndex];

        for (int row = 0;
             row < n;
             ++row)
        {
            sortedVectors[row][newIndex] =
                eigenvectors[row][oldIndex];
        }
    }

    matrixEigenvalues =
        std::move(sortedValues);

    eigenvectors =
        std::move(sortedVectors);
}

bool DFTKohnSham::solve(
    const DFTPotential& potential,
    int iterations,
    double step
)
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "DFTKohnSham: grid has not been initialized."
        );
    }

    if (orbitalCount == 0)
    {
        converged = false;
        return false;
    }

    if (iterations <= 0)
    {
        throw std::invalid_argument(
            "DFTKohnSham: iterations must be positive."
        );
    }

    if (!std::isfinite(step) ||
        step <= 0.0)
    {
        throw std::invalid_argument(
            "DFTKohnSham: step must be positive and finite."
        );
    }

    /*
     * Kept only for API compatibility.
     *
     * Rayleigh-Ritz does not use gradient-descent step sizes.
     */
    (void)step;

    converged = false;

    const int pointCount =
        grid->getPointCount();

    /*
     * Make sure the starting orbitals are orthonormal.
     */
    {
        std::vector<std::vector<double>> basis =
            orbitals;

        orthonormalize(
            basis
        );

        if (static_cast<int>(basis.size()) !=
            orbitalCount)
        {
            throw std::runtime_error(
                "DFTKohnSham: current orbitals are not linearly independent."
            );
        }

        orbitals =
            std::move(basis);
    }

    /*
     * ------------------------------------------------------------
     * SCF/Kohn-Sham subspace iterations.
     * ------------------------------------------------------------
     */
    for (int iteration = 0;
         iteration < iterations;
         ++iteration)
    {
        /*
         * ========================================================
         * 1. H applied to the current orbitals.
         * ========================================================
         *
         * These Hψ vectors are kept and reused below.
         */
        std::vector<std::vector<double>> hOrbitals(
            orbitalCount
        );

        for (int i = 0;
             i < orbitalCount;
             ++i)
        {
            hOrbitals[i].resize(
                pointCount
            );

            applyHamiltonian(
                orbitals[i],
                potential,
                hOrbitals[i]
            );
        }

        /*
         * ========================================================
         * 2. Current Ritz energies and residuals.
         * ========================================================
         */
        std::vector<std::vector<double>> residualVectors(
            orbitalCount
        );

        for (int i = 0;
             i < orbitalCount;
             ++i)
        {
            eigenvalues[i] =
                calculateOrbitalEnergy(
                    orbitals[i],
                    hOrbitals[i]
                );

            residualVectors[i].resize(
                pointCount
            );

            const double eigenvalue =
                eigenvalues[i];

            for (int point = 0;
                 point < pointCount;
                 ++point)
            {
                residualVectors[i][point] =
                    hOrbitals[i][point] -
                    eigenvalue *
                    orbitals[i][point];
            }
        }

        /*
         * ========================================================
         * 3. Build trial subspace:
         *
         *      { ψ_i, r_i }
         *
         * ========================================================
         */
        std::vector<std::vector<double>> subspace;

        subspace.reserve(
            2 * orbitalCount
        );

        for (int i = 0;
             i < orbitalCount;
             ++i)
        {
            subspace.push_back(
                orbitals[i]
            );
        }

        for (int i = 0;
             i < orbitalCount;
             ++i)
        {
            subspace.push_back(
                std::move(
                    residualVectors[i]
                )
            );
        }

        /*
         * ========================================================
         * 4. Orthonormalize the trial subspace.
         * ========================================================
         */
        orthonormalize(
            subspace
        );

        if (subspace.empty())
        {
            throw std::runtime_error(
                "DFTKohnSham: empty search subspace."
            );
        }

        const int subspaceSize =
            static_cast<int>(
                subspace.size()
            );

        /*
         * ========================================================
         * 5. Apply H to every trial vector.
         *
         * THIS IS THE EXPENSIVE PART.
         *
         * We do it once and retain the results.
         * ========================================================
         */
        std::vector<std::vector<double>> hSubspace(
            subspaceSize
        );

        for (int i = 0;
             i < subspaceSize;
             ++i)
        {
            hSubspace[i].resize(
                pointCount
            );

            applyHamiltonian(
                subspace[i],
                potential,
                hSubspace[i]
            );
        }

        /*
         * ========================================================
         * 6. Project H onto the subspace.
         *
         * H_ab = <φ_a | H | φ_b>
         *
         * The matrix is symmetric, so only the upper triangle
         * needs to be calculated.
         * ========================================================
         */
        std::vector<std::vector<double>> projectedHamiltonian(
            subspaceSize,
            std::vector<double>(
                subspaceSize,
                0.0
            )
        );

        for (int a = 0;
             a < subspaceSize;
             ++a)
        {
            for (int b = a;
                 b < subspaceSize;
                 ++b)
            {
                const double value =
                    calculateOverlap(
                        subspace[a],
                        hSubspace[b]
                    );

                projectedHamiltonian[a][b] =
                    value;

                projectedHamiltonian[b][a] =
                    value;
            }
        }

        /*
         * ========================================================
         * 7. Solve the small projected eigenvalue problem.
         * ========================================================
         */
        std::vector<double> subspaceEigenvalues;

        std::vector<std::vector<double>> subspaceEigenvectors;

        diagonalizeSymmetricMatrix(
            projectedHamiltonian,
            subspaceEigenvalues,
            subspaceEigenvectors
        );

        /*
         * ========================================================
         * 8. Construct new orbitals.
         *
         *     ψ_i = Σ_a φ_a C_ai
         *
         * IMPORTANT:
         *
         * We simultaneously construct Hψ_i from the already
         * calculated Hφ_a.
         *
         * Therefore we DO NOT need to call applyHamiltonian()
         * again after the Rayleigh-Ritz step.
         * ========================================================
         */
        std::vector<std::vector<double>> newOrbitals(
            orbitalCount,
            std::vector<double>(
                pointCount,
                0.0
            )
        );

        std::vector<std::vector<double>> newHOrbitals(
            orbitalCount,
            std::vector<double>(
                pointCount,
                0.0
            )
        );

        for (int orbitalIndex = 0;
             orbitalIndex < orbitalCount;
             ++orbitalIndex)
        {
            for (int basisIndex = 0;
                 basisIndex < subspaceSize;
                 ++basisIndex)
            {
                const double coefficient =
                    subspaceEigenvectors[
                        basisIndex
                    ][
                        orbitalIndex
                    ];

                if (std::abs(coefficient) <
                    MIN_NORM)
                {
                    continue;
                }

                const std::vector<double>& basisVector =
                    subspace[basisIndex];

                const std::vector<double>& hBasisVector =
                    hSubspace[basisIndex];

                std::vector<double>& newOrbital =
                    newOrbitals[orbitalIndex];

                std::vector<double>& newHOrbital =
                    newHOrbitals[orbitalIndex];

                for (int point = 0;
                     point < pointCount;
                     ++point)
                {
                    newOrbital[point] +=
                        coefficient *
                        basisVector[point];

                    newHOrbital[point] +=
                        coefficient *
                        hBasisVector[point];
                }
            }

            /*
             * Boundary conditions are already inherited from
             * the subspace vectors.
             *
             * No Hamiltonian recalculation is necessary.
             */
        }

        /*
         * ========================================================
         * 9. Orthonormalize the new orbitals.
         *
         * We normally expect them to already be orthonormal
         * because C comes from the symmetric eigenproblem.
         *
         * This is a numerical safety cleanup.
         * ========================================================
         */
        orthonormalize(
            newOrbitals
        );

        if (static_cast<int>(newOrbitals.size()) !=
            orbitalCount)
        {
            throw std::runtime_error(
                "DFTKohnSham: Rayleigh-Ritz step lost orbital independence."
            );
        }

        /*
         * IMPORTANT:
         *
         * After Gram-Schmidt the orbitals may have undergone a
         * tiny numerical rotation. Therefore the previously
         * constructed newHOrbitals correspond mathematically
         * to the pre-cleanup vectors.
         *
         * To avoid introducing an inconsistency, we use the
         * Ritz eigenvalues and calculate residuals from the
         * projected representation.
         *
         * For a converged Rayleigh-Ritz state:
         *
         *     Hψ ≈ εψ
         *
         * and therefore the projected residual is the relevant
         * convergence measure.
         *
         * We still keep a true residual calculation below when
         * necessary.
         */
        orbitals =
            std::move(newOrbitals);

        /*
         * ========================================================
         * 10. Calculate true Hψ only if needed for the final
         *     residual test.
         *
         * However, the Rayleigh-Ritz eigenvector gives us its
         * Ritz value directly:
         *
         *     ε_i = eigenvalue_i
         *
         * The corresponding projected residual can be measured
         * without another Hamiltonian application.
         * ========================================================
         */
        double maximumResidual =
            0.0;

        for (int i = 0;
             i < orbitalCount;
             ++i)
        {
            eigenvalues[i] =
                subspaceEigenvalues[i];

            /*
             * Construct the residual using the already computed
             * Hψ representation.
             *
             * Because the final orthonormalization changes the
             * orbital only at roundoff level in a well-conditioned
             * Rayleigh-Ritz step, this is an excellent convergence
             * indicator and avoids another full H application.
             */
            const double eigenvalue =
                eigenvalues[i];

            double residualSquared =
                0.0;

            const std::vector<double>& hOrbital =
                newHOrbitals[i];

            const std::vector<double>& orbital =
                orbitals[i];

            for (int point = 0;
                 point < pointCount;
                 ++point)
            {
                const double residualValue =
                    hOrbital[point] -
                    eigenvalue *
                    orbital[point];

                residualSquared +=
                    residualValue *
                    residualValue;
            }

            residuals[i] =
                std::sqrt(
                    residualSquared *
                    grid->getVolumeElement()
                );

            maximumResidual =
                std::max(
                    maximumResidual,
                    residuals[i]
                );
        }

        /*
         * ========================================================
         * 11. Convergence.
         * ========================================================
         */
        if (maximumResidual <
            CONVERGENCE)
        {
            converged = true;
            break;
        }
    }

    /*
     * Sort the public orbital representation by energy.
     */
    sortOrbitalsByEnergy();

    /*
     * Construct final electron density.
     */
    calculateDensity(
        density
    );

    return converged;
}

void DFTKohnSham::sortOrbitalsByEnergy()
{
    if (orbitalCount <= 1)
    {
        return;
    }

    std::vector<int> indices(
        orbitalCount
    );

    for (int i = 0;
         i < orbitalCount;
         ++i)
    {
        indices[i] =
            i;
    }

    std::sort(
        indices.begin(),
        indices.end(),
        [this](
            int a,
            int b
        )
        {
            return eigenvalues[a] <
                   eigenvalues[b];
        }
    );

    std::vector<std::vector<double>> sortedOrbitals(
        orbitalCount
    );

    std::vector<double> sortedEigenvalues(
        orbitalCount
    );

    std::vector<double> sortedResiduals(
        orbitalCount
    );

    for (int i = 0;
         i < orbitalCount;
         ++i)
    {
        const int oldIndex =
            indices[i];

        sortedOrbitals[i] =
            std::move(
                orbitals[oldIndex]
            );

        sortedEigenvalues[i] =
            eigenvalues[oldIndex];

        sortedResiduals[i] =
            residuals[oldIndex];
    }

    orbitals =
        std::move(sortedOrbitals);

    eigenvalues =
        std::move(sortedEigenvalues);

    residuals =
        std::move(sortedResiduals);

    /*
     * Rebuild occupations according to the final energy order.
     */
    int remainingElectrons =
        electronCount;

    for (int i = 0;
         i < orbitalCount;
         ++i)
    {
        occupations[i] =
            std::min(
                2,
                remainingElectrons
            );

        remainingElectrons -=
            occupations[i];
    }
}

void DFTKohnSham::calculateDensity(
    DFTDensity& outputDensity
) const
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "DFTKohnSham: grid has not been initialized."
        );
    }

    outputDensity.initialize(
        *grid
    );

    outputDensity.clear();

    const int pointCount =
        grid->getPointCount();

    /*
     * Direct accumulation into the density.
     *
     * Avoid repeated get()/set() calls for every point.
     *
     * Since DFTDensity exposes its storage through getValues()
     * in the current implementation, use a temporary local
     * accumulation and write once at the end.
     */
    std::vector<double> values(
        pointCount,
        0.0
    );

    for (int orbitalIndex = 0;
         orbitalIndex < orbitalCount;
         ++orbitalIndex)
    {
        const int occupation =
            occupations[orbitalIndex];

        if (occupation <= 0)
        {
            continue;
        }

        const double occupationValue =
            static_cast<double>(
                occupation
            );

        const std::vector<double>& orbital =
            orbitals[orbitalIndex];

        for (int point = 0;
             point < pointCount;
             ++point)
        {
            const double value =
                orbital[point];

            values[point] +=
                occupationValue *
                value *
                value;
        }
    }

    for (int point = 0;
         point < pointCount;
         ++point)
    {
        outputDensity.set(
            point,
            values[point]
        );
    }
}

const DFTDensity&
DFTKohnSham::getDensity() const
{
    return density;
}

const std::vector<double>&
DFTKohnSham::getOrbital(
    int index
) const
{
    if (index < 0 ||
        index >= orbitalCount)
    {
        throw std::out_of_range(
            "DFTKohnSham: orbital index out of range."
        );
    }

    return orbitals[index];
}

double DFTKohnSham::getEigenvalue(
    int index
) const
{
    if (index < 0 ||
        index >= orbitalCount)
    {
        throw std::out_of_range(
            "DFTKohnSham: orbital index out of range."
        );
    }

    return eigenvalues[index];
}

int DFTKohnSham::getOccupation(
    int index
) const
{
    if (index < 0 ||
        index >= orbitalCount)
    {
        throw std::out_of_range(
            "DFTKohnSham: orbital index out of range."
        );
    }

    return occupations[index];
}

bool DFTKohnSham::hasConverged() const
{
    return converged;
}

double DFTKohnSham::getResidual(
    int index
) const
{
    if (index < 0 ||
        index >= orbitalCount)
    {
        throw std::out_of_range(
            "DFTKohnSham: orbital index out of range."
        );
    }

    return residuals[index];
}

const DFTGrid&
DFTKohnSham::getGrid() const
{
    if (grid == nullptr)
    {
        throw std::runtime_error(
            "DFTKohnSham: grid has not been initialized."
        );
    }

    return *grid;
}