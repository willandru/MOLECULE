#pragma once


namespace GeometryConstants
{
    // ========================================================
    // SPATIAL UNIT
    // ========================================================

    // 1 internal geometry unit = 1 Angstrom.
    constexpr float ANGSTROM = 1.0f;

    constexpr float NANOMETER =
        10.0f * ANGSTROM;

    constexpr float MICROMETER =
        10000.0f * ANGSTROM;


    // ========================================================
    // ATOMIC UNITS
    // ========================================================

    // 1 Bohr radius in Angstrom.
    //
    // a0 = 0.529177210903 Å

    constexpr double BOHR_TO_ANGSTROM =
        0.529177210903;

    constexpr double ANGSTROM_TO_BOHR =
        1.0 / BOHR_TO_ANGSTROM;


    // ========================================================
    // GRID
    // ========================================================

    // Reference grid size:
    //
    // 100 Å × 100 Å × 100 Å

    constexpr float GRID_SIZE =
        100.0f * ANGSTROM;

    constexpr int GRID_DIVISIONS =
        100;


    // ========================================================
    // GRID SPACING
    // ========================================================

    constexpr float GRID_SPACING =
        GRID_SIZE /
        static_cast<float>(
            GRID_DIVISIONS
        );
}