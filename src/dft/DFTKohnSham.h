#pragma once

#include <vector>

#include "DFTGrid.h"
#include "DFTPotential.h"
#include "DFTDensity.h"

class DFTKohnSham
{
public:
    DFTKohnSham();
    explicit DFTKohnSham(const DFTGrid& grid);

    void initialize(const DFTGrid& grid);

    void setElectronCount(int count);
    int getElectronCount() const;

    int getOrbitalCount() const;

    bool solve(
        const DFTPotential& potential,
        int iterations = 200,
        double step = 0.001
    );

    void calculateDensity(DFTDensity& density) const;

    const DFTDensity& getDensity() const;

    const std::vector<double>& getOrbital(int index) const;

    double getEigenvalue(int index) const;
    int getOccupation(int index) const;

    bool hasConverged() const;

    double getResidual(int index) const;

    const DFTGrid& getGrid() const;

private:
    const DFTGrid* grid;

    int electronCount;
    int orbitalCount;

    std::vector<std::vector<double>> orbitals;
    std::vector<double> eigenvalues;
    std::vector<int> occupations;
    std::vector<double> residuals;

    DFTDensity density;

    bool converged;

    void initializeOrbitals();

    double calculateNorm(
        const std::vector<double>& orbital
    ) const;

    double calculateOverlap(
        const std::vector<double>& a,
        const std::vector<double>& b
    ) const;

    void orthonormalize(
        std::vector<std::vector<double>>& basis
    ) const;

    void applyOrbitalBoundaryConditions(
        std::vector<double>& orbital
    ) const;

    void calculateLaplacian(
        const std::vector<double>& orbital,
        std::vector<double>& laplacian
    ) const;

    void applyHamiltonian(
        const std::vector<double>& orbital,
        const DFTPotential& potential,
        std::vector<double>& result
    ) const;

    double calculateOrbitalEnergy(
        const std::vector<double>& orbital,
        const std::vector<double>& hOrbital
    ) const;

    double calculateResidual(
        const std::vector<double>& orbital,
        const std::vector<double>& hOrbital,
        double eigenvalue
    ) const;

    void diagonalizeSymmetricMatrix(
        std::vector<std::vector<double>>& matrix,
        std::vector<double>& eigenvalues,
        std::vector<std::vector<double>>& eigenvectors
    ) const;

    void sortOrbitalsByEnergy();
};