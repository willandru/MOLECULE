#include "DFT.h"
#include "DFTGrid.h"

#include <glm/glm.hpp>

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>


// ================================================================
// UTILIDADES
// ================================================================

void printEnergyReport(
    const std::string& name,
    const DFT& dft
)
{
    std::cout << "\n";
    std::cout
        << "============================================================\n";

    std::cout
        << name
        << "\n";

    std::cout
        << "============================================================\n";

    std::cout
        << "SCF converged: "
        << (dft.hasConverged() ? "YES" : "NO")
        << "\n";

    std::cout
        << "Electrons: "
        << dft.getKohnSham().getElectronCount()
        << "\n";

    std::cout
        << "Orbitals: "
        << dft.getKohnSham().getOrbitalCount()
        << "\n";

    std::cout
        << std::setprecision(12);

    std::cout
        << "E electron-nucleus : "
        << dft.getElectronNuclearEnergy()
        << " Ha\n";

    std::cout
        << "E Hartree          : "
        << dft.getHartreeEnergy()
        << " Ha\n";

    std::cout
        << "E exchange         : "
        << dft.getExchangeEnergy()
        << " Ha\n";

    std::cout
        << "E nucleus-nucleus  : "
        << dft.getNuclearRepulsionEnergy()
        << " Ha\n";

    std::cout
        << "E TOTAL            : "
        << dft.calculateTotalEnergy()
        << " Ha\n";

    std::cout
        << "\nKohn-Sham orbitals:\n";

    for (
        int i = 0;
        i < dft.getKohnSham().getOrbitalCount();
        ++i
    )
    {
        std::cout
            << "  Orbital "
            << i
            << " | occupation = "
            << dft.getKohnSham().getOccupation(i)
            << " | eigenvalue = "
            << dft.getKohnSham().getEigenvalue(i)
            << " Ha"
            << " | residual = "
            << dft.getKohnSham().getResidual(i)
            << "\n";
    }
}


// ================================================================
// ATOM TEST
// ================================================================

void runAtomTest(
    const std::string& name,
    int atomicNumber,
    const DFTGrid& grid
)
{
    std::cout << "\n\n";

    std::cout
        << "############################################################\n";

    std::cout
        << "# "
        << name
        << "\n";

    std::cout
        << "############################################################\n";

    std::cout
        << "["
        << name
        << "] Creando objeto DFT..."
        << std::endl;

    DFT atom(grid);

    std::cout
        << "["
        << name
        << "] Objeto DFT creado."
        << std::endl;


    // ------------------------------------------------------------
    // NÚCLEO
    // ------------------------------------------------------------

    std::vector<int> charges =
    {
        atomicNumber
    };

    std::vector<glm::dvec3> positions =
    {
        glm::dvec3(0.0)
    };

    std::cout
        << "["
        << name
        << "] Antes de setMolecule()"
        << std::endl;

    atom.setMolecule(
        charges,
        positions
    );

    std::cout
        << "["
        << name
        << "] Despues de setMolecule()"
        << std::endl;


    // ------------------------------------------------------------
    // POTENCIAL
    // ------------------------------------------------------------

    atom.setPotentialSoftening(
        0.15
    );

    std::cout
        << "["
        << name
        << "] Softening configurado."
        << std::endl;


    // ------------------------------------------------------------
    // SCF
    // ------------------------------------------------------------

    std::cout
        << "["
        << name
        << "] Antes de solveSCF()"
        << std::endl;

    const bool converged =
        atom.solveSCF(
            60,
            200,
            0.001,
            0.20,
            1.0e-5
        );

    std::cout
        << "["
        << name
        << "] Despues de solveSCF(): "
        << (
            converged
            ? "CONVERGED"
            : "NO CONVERGED"
        )
        << std::endl;


    // ------------------------------------------------------------
    // RESULTADOS
    // ------------------------------------------------------------

    printEnergyReport(
        name,
        atom
    );
}


// ================================================================
// MAIN
// ================================================================

int main()
{
    std::cout
        << "\n";

    std::cout
        << "============================================================\n";

    std::cout
        << "DFT ATOMIC SIMULATION\n";

    std::cout
        << "============================================================\n";


    // ============================================================
    // DFT GRID
    // ============================================================

    /*
        Unidades atómicas:

            distancia = Bohr
            energía   = Hartree

        Malla espacial:

            32 x 32 x 32 puntos
            spacing = 0.5 Bohr
            box = 16 Bohr
    */

    const int NX = 32;
    const int NY = 32;
    const int NZ = 32;

    const double DFT_SPACING = 0.5;

    const double DFT_BOX =
        DFT_SPACING *
        static_cast<double>(NX);

    const glm::dvec3 DFT_ORIGIN(
        -0.5 * DFT_BOX,
        -0.5 * DFT_BOX,
        -0.5 * DFT_BOX
    );

    DFTGrid dftGrid(
        NX,
        NY,
        NZ,
        DFT_SPACING,
        DFT_ORIGIN
    );


    std::cout
        << "Grid: "
        << NX
        << " x "
        << NY
        << " x "
        << NZ
        << "\n";

    std::cout
        << "Spacing: "
        << DFT_SPACING
        << " Bohr\n";

    std::cout
        << "Box size: "
        << DFT_BOX
        << " Bohr\n";


    // ============================================================
    // H
    // ============================================================

    runAtomTest(
        "H",
        1,
        dftGrid
    );


    // ============================================================
    // He
    // ============================================================

    runAtomTest(
        "He",
        2,
        dftGrid
    );


    // ============================================================
    // Li
    // ============================================================

    runAtomTest(
        "Li",
        3,
        dftGrid
    );


    // ============================================================
    // FIN
    // ============================================================

    std::cout
        << "\n";

    std::cout
        << "============================================================\n";

    std::cout
        << "DFT ATOMIC TESTS FINISHED\n";

    std::cout
        << "============================================================\n";


    return 0;
}