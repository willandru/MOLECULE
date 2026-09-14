#include "Window.h"
#include "Camera.h"
#include "InputKeyboard.h"
#include "InputMouse.h"
#include "Timer1.h"

#include "Grid.h"
#include "GridRenderer.h"

#include "Atom.h"
#include "AtomRenderer.h"

#include <glad/glad.h>

#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <array>


int main()
{
    Window window(
        1280,
        720,
        "PERIODIC TABLE",
        true
    );


    // ========================================================
    // CAMERA
    // ========================================================

    Camera camera(
        glm::vec3(0.0f),
        70.0f
    );

    camera.setAspectRatio(
        window.getAspectRatio()
    );

    camera.setOrientation(
        0.0f,
        glm::radians(-70.0f)
    );


    // ========================================================
    // INPUT
    // ========================================================

    InputKeyboard keyboard(
        window.getHandle()
    );

    InputMouse mouse(
        window.getHandle()
    );


    // ========================================================
    // TIMER
    // ========================================================

    Timer1 timer;


    // ========================================================
    // GRID
    // ========================================================

    Grid grid;

    GridRenderer gridRenderer;

    gridRenderer.initialize(
        grid
    );


    // ========================================================
    // PERIODIC TABLE
    // ========================================================
    //
    // 18 GROUPS
    //
    //                    GROUP
    //
    //       1   2   3 ... 12  13 14 15 16 17 18
    //
    // Period 1:
    //       H                               He
    //
    // Period 2:
    //       Li  Be                  B  C  N  O  F  Ne
    //
    // Period 3:
    //       Na  Mg                 Al Si P  S  Cl Ar
    //
    // Period 4:
    //       K   Ca  Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr
    //
    // Period 5:
    //       Rb  Sr  Y  Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe
    //
    // Period 6:
    //       Cs  Ba  La Hf Ta W  Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn
    //
    // Period 7:
    //       Fr  Ra  Ac Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og
    //
    //
    // LANTHANIDES:
    //
    //       Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu
    //
    // ACTINIDES:
    //
    //       Th Pa U  Np Pu Am Cm Bk Cf Es Fm Md No Lr
    //
    // ========================================================

    constexpr int ROWS = 7;
    constexpr int GROUPS = 18;

    constexpr int periodicTable[ROWS][GROUPS] =
    {
        // ====================================================
        // PERIOD 1
        // ====================================================

        {
            1,  0,  0,  0,  0,  0,  0,  0,  0,
            0,  0,  0,  0,  0,  0,  0,  0,  2
        },


        // ====================================================
        // PERIOD 2
        // ====================================================

        {
            3,  4,  0,  0,  0,  0,  0,  0,  0,
            0,  0,  0,  5,  6,  7,  8,  9, 10
        },


        // ====================================================
        // PERIOD 3
        // ====================================================

        {
            11, 12, 0,  0,  0,  0,  0,  0,  0,
            0,  0,  0, 13, 14, 15, 16, 17, 18
        },


        // ====================================================
        // PERIOD 4
        // ====================================================

        {
            19, 20, 21, 22, 23, 24, 25, 26, 27,
            28, 29, 30, 31, 32, 33, 34, 35, 36
        },


        // ====================================================
        // PERIOD 5
        // ====================================================

        {
            37, 38, 39, 40, 41, 42, 43, 44, 45,
            46, 47, 48, 49, 50, 51, 52, 53, 54
        },


        // ====================================================
        // PERIOD 6
        // ====================================================

        {
            55, 56, 57, 72, 73, 74, 75, 76, 77,
            78, 79, 80, 81, 82, 83, 84, 85, 86
        },


        // ====================================================
        // PERIOD 7
        // ====================================================

        {
            87, 88, 89, 104, 105, 106, 107, 108, 109,
            110, 111, 112, 113, 114, 115, 116, 117, 118
        }
    };


    // ========================================================
    // F-BLOCK
    // ========================================================
    //
    // Lanthanides:
    //
    // Ce - Lu = 58 - 71
    //
    // Actinides:
    //
    // Th - Lr = 90 - 103
    //
    // ========================================================

    constexpr int lanthanides[14] =
    {
        58, 59, 60, 61, 62, 63, 64,
        65, 66, 67, 68, 69, 70, 71
    };


    constexpr int actinides[14] =
    {
        90, 91, 92, 93, 94, 95, 96,
        97, 98, 99, 100, 101, 102, 103
    };


    // ========================================================
    // GEOMETRY
    // ========================================================

    constexpr float spacingX = 3.5f;
    constexpr float spacingY = 3.5f;

    constexpr float fBlockOffsetX = 3.5f;

    constexpr float tableWidth =
        static_cast<float>(GROUPS - 1) * spacingX;

    constexpr float tableCenterX =
        tableWidth * 0.5f;


    // ========================================================
    // CREATE ATOMS
    // ========================================================

    std::vector<Atom> atoms;

    atoms.reserve(118);


    // ========================================================
    // MAIN PERIODIC TABLE
    // ========================================================

    for (int period = 0; period < ROWS; ++period)
    {
        for (int group = 0; group < GROUPS; ++group)
        {
            const int atomicNumber =
                periodicTable[period][group];


            if (atomicNumber == 0)
            {
                continue;
            }


            const float x =
                static_cast<float>(group) *
                spacingX -
                tableCenterX;


            const float y =
                -static_cast<float>(period) *
                spacingY +
                10.5f;


            atoms.emplace_back(
                atomicNumber,
                glm::vec3(
                    x,
                    y,
                    0.0f
                )
            );
        }
    }


    // ========================================================
    // LANTHANIDES
    // ========================================================

    for (int i = 0; i < 14; ++i)
    {
        const float x =
            static_cast<float>(i + 2) *
            spacingX -
            tableCenterX;


        const float y =
            10.5f -
            7.0f * spacingY;


        atoms.emplace_back(
            lanthanides[i],
            glm::vec3(
                x,
                y,
                0.0f
            )
        );
    }


    // ========================================================
    // ACTINIDES
    // ========================================================

    for (int i = 0; i < 14; ++i)
    {
        const float x =
            static_cast<float>(i + 2) *
            spacingX -
            tableCenterX;


        const float y =
            10.5f -
            8.0f * spacingY;


        atoms.emplace_back(
            actinides[i],
            glm::vec3(
                x,
                y,
                0.0f
            )
        );
    }


    // ========================================================
    // ATOM RENDERER
    // ========================================================

    AtomRenderer atomRenderer;

    if (!atomRenderer.initialize())
    {
        return -1;
    }


    // ========================================================
    // OPENGL
    // ========================================================

    glEnable(GL_DEPTH_TEST);


    // ========================================================
    // MAIN LOOP
    // ========================================================

    while (!window.shouldClose())
    {
        timer.update();

        window.pollEvents();


        // ----------------------------------------------------
        // CLOSE
        // ----------------------------------------------------

        if (keyboard.shouldClose())
        {
            break;
        }


        // ----------------------------------------------------
        // KEYBOARD
        // ----------------------------------------------------

        keyboard.update(
            camera,
            timer.getDeltaTime()
        );


        // ----------------------------------------------------
        // MOUSE
        // ----------------------------------------------------

        mouse.update();


        if (mouse.isMiddleButtonPressed())
        {
            const glm::vec2& delta =
                mouse.getDelta();

            camera.orbit(
                delta.x,
                delta.y
            );
        }


        // ----------------------------------------------------
        // ZOOM
        // ----------------------------------------------------

        const float scroll =
            mouse.getScrollDelta();


        if (scroll != 0.0f)
        {
            camera.zoom(
                scroll
            );
        }


        mouse.clearScrollDelta();


        // ====================================================
        // CLEAR
        // ====================================================

        glClear(
            GL_COLOR_BUFFER_BIT |
            GL_DEPTH_BUFFER_BIT
        );


        // ====================================================
        // GRID
        // ====================================================

        gridRenderer.render(
            camera.getViewMatrix(),
            camera.getProjectionMatrix()
        );


        // ====================================================
        // ATOMS
        // ====================================================

        for (const Atom& atom : atoms)
        {
            atomRenderer.render(
                atom,
                camera.getViewMatrix(),
                camera.getProjectionMatrix()
            );
        }


        // ====================================================
        // PRESENT
        // ====================================================

        window.swapBuffers();
    }


    // ========================================================
    // CLEANUP
    // ========================================================

    atomRenderer.shutdown();


    return 0;
}