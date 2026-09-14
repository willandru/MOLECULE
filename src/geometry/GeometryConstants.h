#pragma once


namespace GeometryConstants
{
    // ========================================================
    // SPATIAL UNIT
    // ========================================================

    // 1 internal geometry unit = 1 Angstrom.
    constexpr float ANGSTROM = 1.0f;

    constexpr float NANOMETER = 10.0f * ANGSTROM;
    constexpr float MICROMETER = 10000.0f * ANGSTROM;


    // ========================================================
    // GRID
    // ========================================================

    // Reference grid size:
    //
    // 10 Å × 10 Å × 10 Å
    //
    // This gives a convenient visual reference for
    // atomic and small molecular structures.

    constexpr float GRID_SIZE = 10.0f * ANGSTROM;

    constexpr int GRID_DIVISIONS = 10;


    // ========================================================
    // GRID SPACING
    // ========================================================

    constexpr float GRID_SPACING =
        GRID_SIZE /
        static_cast<float>(GRID_DIVISIONS);
}