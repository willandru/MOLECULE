#pragma once

#include "GeometryConstants.h"

// ============================================================
// COVALENT RADIUS CONSTANTS
// ============================================================
//
// Covalent radii for the basic atomic representation.
//
// Reference:
//     P. Pyykkö and M. Atsumi,
//     "Molecular Single-Bond Covalent Radii for Elements 1–118",
//     Chemistry – A European Journal, 15, 186–197 (2009).
//
// The values correspond to single-bond covalent radii.
//
// Source values are defined in picometers (pm).
//
// Internal geometry unit:
//     1.0 = 1 Angstrom (Å)
//
// Therefore:
//
//     1 Å = 100 pm
//
// The values returned by get() are converted to the internal
// geometry unit (Å).
//
// ============================================================

namespace CovalentRadius
{
    // ========================================================
    // CONVERSION
    // ========================================================

    constexpr float PM_TO_ANGSTROM = 0.01f;


    // ========================================================
    // PERIOD 1
    // ========================================================

    constexpr float H  = 32.0f * PM_TO_ANGSTROM;
    constexpr float He = 46.0f * PM_TO_ANGSTROM;


    // ========================================================
    // PERIOD 2
    // ========================================================

    constexpr float Li = 133.0f * PM_TO_ANGSTROM;
    constexpr float Be = 102.0f * PM_TO_ANGSTROM;
    constexpr float B  = 85.0f * PM_TO_ANGSTROM;
    constexpr float C  = 75.0f * PM_TO_ANGSTROM;
    constexpr float N  = 71.0f * PM_TO_ANGSTROM;
    constexpr float O  = 63.0f * PM_TO_ANGSTROM;
    constexpr float F  = 64.0f * PM_TO_ANGSTROM;
    constexpr float Ne = 67.0f * PM_TO_ANGSTROM;


    // ========================================================
    // PERIOD 3
    // ========================================================

    constexpr float Na = 155.0f * PM_TO_ANGSTROM;
    constexpr float Mg = 139.0f * PM_TO_ANGSTROM;
    constexpr float Al = 126.0f * PM_TO_ANGSTROM;
    constexpr float Si = 116.0f * PM_TO_ANGSTROM;
    constexpr float P  = 111.0f * PM_TO_ANGSTROM;
    constexpr float S  = 103.0f * PM_TO_ANGSTROM;
    constexpr float Cl = 99.0f  * PM_TO_ANGSTROM;
    constexpr float Ar = 96.0f  * PM_TO_ANGSTROM;


    // ========================================================
    // PERIOD 4
    // ========================================================

    constexpr float K  = 196.0f * PM_TO_ANGSTROM;
    constexpr float Ca = 171.0f * PM_TO_ANGSTROM;
    constexpr float Sc = 148.0f * PM_TO_ANGSTROM;
    constexpr float Ti = 136.0f * PM_TO_ANGSTROM;
    constexpr float V  = 134.0f * PM_TO_ANGSTROM;
    constexpr float Cr = 122.0f * PM_TO_ANGSTROM;
    constexpr float Mn = 119.0f * PM_TO_ANGSTROM;
    constexpr float Fe = 116.0f * PM_TO_ANGSTROM;
    constexpr float Co = 111.0f * PM_TO_ANGSTROM;
    constexpr float Ni = 110.0f * PM_TO_ANGSTROM;
    constexpr float Cu = 112.0f * PM_TO_ANGSTROM;
    constexpr float Zn = 118.0f * PM_TO_ANGSTROM;
    constexpr float Ga = 124.0f * PM_TO_ANGSTROM;
    constexpr float Ge = 121.0f * PM_TO_ANGSTROM;
    constexpr float As = 121.0f * PM_TO_ANGSTROM;
    constexpr float Se = 116.0f * PM_TO_ANGSTROM;
    constexpr float Br = 114.0f * PM_TO_ANGSTROM;
    constexpr float Kr = 117.0f * PM_TO_ANGSTROM;


    // ========================================================
    // PERIOD 5
    // ========================================================

    constexpr float Rb = 210.0f * PM_TO_ANGSTROM;
    constexpr float Sr = 185.0f * PM_TO_ANGSTROM;
    constexpr float Y  = 163.0f * PM_TO_ANGSTROM;
    constexpr float Zr = 154.0f * PM_TO_ANGSTROM;
    constexpr float Nb = 147.0f * PM_TO_ANGSTROM;
    constexpr float Mo = 138.0f * PM_TO_ANGSTROM;
    constexpr float Tc = 128.0f * PM_TO_ANGSTROM;
    constexpr float Ru = 125.0f * PM_TO_ANGSTROM;
    constexpr float Rh = 125.0f * PM_TO_ANGSTROM;
    constexpr float Pd = 120.0f * PM_TO_ANGSTROM;
    constexpr float Ag = 128.0f * PM_TO_ANGSTROM;
    constexpr float Cd = 136.0f * PM_TO_ANGSTROM;
    constexpr float In = 142.0f * PM_TO_ANGSTROM;
    constexpr float Sn = 140.0f * PM_TO_ANGSTROM;
    constexpr float Sb = 140.0f * PM_TO_ANGSTROM;
    constexpr float Te = 136.0f * PM_TO_ANGSTROM;
    constexpr float I  = 133.0f * PM_TO_ANGSTROM;
    constexpr float Xe = 131.0f * PM_TO_ANGSTROM;


    // ========================================================
    // PERIOD 6
    // ========================================================

    constexpr float Cs = 232.0f * PM_TO_ANGSTROM;
    constexpr float Ba = 196.0f * PM_TO_ANGSTROM;

    constexpr float La = 180.0f * PM_TO_ANGSTROM;
    constexpr float Ce = 163.0f * PM_TO_ANGSTROM;
    constexpr float Pr = 176.0f * PM_TO_ANGSTROM;
    constexpr float Nd = 174.0f * PM_TO_ANGSTROM;
    constexpr float Pm = 173.0f * PM_TO_ANGSTROM;
    constexpr float Sm = 172.0f * PM_TO_ANGSTROM;
    constexpr float Eu = 168.0f * PM_TO_ANGSTROM;
    constexpr float Gd = 169.0f * PM_TO_ANGSTROM;
    constexpr float Tb = 168.0f * PM_TO_ANGSTROM;
    constexpr float Dy = 167.0f * PM_TO_ANGSTROM;
    constexpr float Ho = 166.0f * PM_TO_ANGSTROM;
    constexpr float Er = 165.0f * PM_TO_ANGSTROM;
    constexpr float Tm = 164.0f * PM_TO_ANGSTROM;
    constexpr float Yb = 170.0f * PM_TO_ANGSTROM;
    constexpr float Lu = 162.0f * PM_TO_ANGSTROM;

    constexpr float Hf = 152.0f * PM_TO_ANGSTROM;
    constexpr float Ta = 146.0f * PM_TO_ANGSTROM;
    constexpr float W  = 137.0f * PM_TO_ANGSTROM;
    constexpr float Re = 131.0f * PM_TO_ANGSTROM;
    constexpr float Os = 129.0f * PM_TO_ANGSTROM;
    constexpr float Ir = 122.0f * PM_TO_ANGSTROM;
    constexpr float Pt = 123.0f * PM_TO_ANGSTROM;
    constexpr float Au = 124.0f * PM_TO_ANGSTROM;
    constexpr float Hg = 133.0f * PM_TO_ANGSTROM;
    constexpr float Tl = 144.0f * PM_TO_ANGSTROM;
    constexpr float Pb = 144.0f * PM_TO_ANGSTROM;
    constexpr float Bi = 151.0f * PM_TO_ANGSTROM;
    constexpr float Po = 145.0f * PM_TO_ANGSTROM;
    constexpr float At = 147.0f * PM_TO_ANGSTROM;
    constexpr float Rn = 142.0f * PM_TO_ANGSTROM;


    // ========================================================
    // PERIOD 7
    // ========================================================

    constexpr float Fr = 223.0f * PM_TO_ANGSTROM;
    constexpr float Ra = 201.0f * PM_TO_ANGSTROM;

    constexpr float Ac = 186.0f * PM_TO_ANGSTROM;
    constexpr float Th = 175.0f * PM_TO_ANGSTROM;
    constexpr float Pa = 169.0f * PM_TO_ANGSTROM;
    constexpr float U  = 170.0f * PM_TO_ANGSTROM;
    constexpr float Np = 171.0f * PM_TO_ANGSTROM;
    constexpr float Pu = 172.0f * PM_TO_ANGSTROM;
    constexpr float Am = 166.0f * PM_TO_ANGSTROM;
    constexpr float Cm = 166.0f * PM_TO_ANGSTROM;
    constexpr float Bk = 168.0f * PM_TO_ANGSTROM;
    constexpr float Cf = 168.0f * PM_TO_ANGSTROM;
    constexpr float Es = 165.0f * PM_TO_ANGSTROM;
    constexpr float Fm = 167.0f * PM_TO_ANGSTROM;
    constexpr float Md = 173.0f * PM_TO_ANGSTROM;
    constexpr float No = 176.0f * PM_TO_ANGSTROM;
    constexpr float Lr = 161.0f * PM_TO_ANGSTROM;

    constexpr float Rf = 157.0f * PM_TO_ANGSTROM;
    constexpr float Db = 149.0f * PM_TO_ANGSTROM;
    constexpr float Sg = 143.0f * PM_TO_ANGSTROM;
    constexpr float Bh = 141.0f * PM_TO_ANGSTROM;
    constexpr float Hs = 134.0f * PM_TO_ANGSTROM;
    constexpr float Mt = 129.0f * PM_TO_ANGSTROM;
    constexpr float Ds = 128.0f * PM_TO_ANGSTROM;
    constexpr float Rg = 121.0f * PM_TO_ANGSTROM;
    constexpr float Cn = 122.0f * PM_TO_ANGSTROM;
    constexpr float Nh = 136.0f * PM_TO_ANGSTROM;
    constexpr float Fl = 143.0f * PM_TO_ANGSTROM;
    constexpr float Mc = 162.0f * PM_TO_ANGSTROM;
    constexpr float Lv = 175.0f * PM_TO_ANGSTROM;
    constexpr float Ts = 165.0f * PM_TO_ANGSTROM;
    constexpr float Og = 157.0f * PM_TO_ANGSTROM;


    // ========================================================
    // GET COVALENT RADIUS
    // ========================================================

    constexpr float get(int atomicNumber)
    {
        switch (atomicNumber)
        {
            // Period 1
            case 1:   return H;
            case 2:   return He;

            // Period 2
            case 3:   return Li;
            case 4:   return Be;
            case 5:   return B;
            case 6:   return C;
            case 7:   return N;
            case 8:   return O;
            case 9:   return F;
            case 10:  return Ne;

            // Period 3
            case 11:  return Na;
            case 12:  return Mg;
            case 13:  return Al;
            case 14:  return Si;
            case 15:  return P;
            case 16:  return S;
            case 17:  return Cl;
            case 18:  return Ar;

            // Period 4
            case 19:  return K;
            case 20:  return Ca;
            case 21:  return Sc;
            case 22:  return Ti;
            case 23:  return V;
            case 24:  return Cr;
            case 25:  return Mn;
            case 26:  return Fe;
            case 27:  return Co;
            case 28:  return Ni;
            case 29:  return Cu;
            case 30:  return Zn;
            case 31:  return Ga;
            case 32:  return Ge;
            case 33:  return As;
            case 34:  return Se;
            case 35:  return Br;
            case 36:  return Kr;

            // Period 5
            case 37:  return Rb;
            case 38:  return Sr;
            case 39:  return Y;
            case 40:  return Zr;
            case 41:  return Nb;
            case 42:  return Mo;
            case 43:  return Tc;
            case 44:  return Ru;
            case 45:  return Rh;
            case 46:  return Pd;
            case 47:  return Ag;
            case 48:  return Cd;
            case 49:  return In;
            case 50:  return Sn;
            case 51:  return Sb;
            case 52:  return Te;
            case 53:  return I;
            case 54:  return Xe;

            // Period 6
            case 55:  return Cs;
            case 56:  return Ba;
            case 57:  return La;
            case 58:  return Ce;
            case 59:  return Pr;
            case 60:  return Nd;
            case 61:  return Pm;
            case 62:  return Sm;
            case 63:  return Eu;
            case 64:  return Gd;
            case 65:  return Tb;
            case 66:  return Dy;
            case 67:  return Ho;
            case 68:  return Er;
            case 69:  return Tm;
            case 70:  return Yb;
            case 71:  return Lu;
            case 72:  return Hf;
            case 73:  return Ta;
            case 74:  return W;
            case 75:  return Re;
            case 76:  return Os;
            case 77:  return Ir;
            case 78:  return Pt;
            case 79:  return Au;
            case 80:  return Hg;
            case 81:  return Tl;
            case 82:  return Pb;
            case 83:  return Bi;
            case 84:  return Po;
            case 85:  return At;
            case 86:  return Rn;

            // Period 7
            case 87:  return Fr;
            case 88:  return Ra;
            case 89:  return Ac;
            case 90:  return Th;
            case 91:  return Pa;
            case 92:  return U;
            case 93:  return Np;
            case 94:  return Pu;
            case 95:  return Am;
            case 96:  return Cm;
            case 97:  return Bk;
            case 98:  return Cf;
            case 99:  return Es;
            case 100: return Fm;
            case 101: return Md;
            case 102: return No;
            case 103: return Lr;
            case 104: return Rf;
            case 105: return Db;
            case 106: return Sg;
            case 107: return Bh;
            case 108: return Hs;
            case 109: return Mt;
            case 110: return Ds;
            case 111: return Rg;
            case 112: return Cn;
            case 113: return Nh;
            case 114: return Fl;
            case 115: return Mc;
            case 116: return Lv;
            case 117: return Ts;
            case 118: return Og;

            default:
                return 0.0f;
        }
    }
}