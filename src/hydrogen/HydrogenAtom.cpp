#include "HydrogenAtom.h"


HydrogenAtom::HydrogenAtom()
    : position(0.0f),
      dft(),
      dftCalculated(false)
{
}


HydrogenAtom::HydrogenAtom(
    const glm::vec3& position
)
    : position(position),
      dft(),
      dftCalculated(false)
{
}


// ============================================================
// POSITION
// ============================================================

void HydrogenAtom::setPosition(
    const glm::vec3& position
)
{
    this->position = position;
}


const glm::vec3& HydrogenAtom::getPosition() const
{
    return position;
}


// ============================================================
// ATOMIC DATA
// ============================================================

int HydrogenAtom::getAtomicNumber() const
{
    return ATOMIC_NUMBER;
}


int HydrogenAtom::getElectronCount() const
{
    return ELECTRON_COUNT;
}


// ============================================================
// DFT
// ============================================================

void HydrogenAtom::calculateDFT()
{
    dft.solve();

    dftCalculated = true;
}


bool HydrogenAtom::isDFTConverged() const
{
    if (!dftCalculated)
    {
        return false;
    }

    return dft.getResult().converged;
}


std::size_t HydrogenAtom::getSCFIterations() const
{
    if (!dftCalculated)
    {
        return 0;
    }

    return dft.getResult().iterations;
}


// ============================================================
// ENERGY
// ============================================================

double HydrogenAtom::getEnergy() const
{
    return dft.getResult().totalEnergy;
}


double HydrogenAtom::getKSEigenvalue() const
{
    return dft.getResult().eigenvalue;
}


double HydrogenAtom::getKineticEnergy() const
{
    return dft.getResult().kineticEnergy;
}


double HydrogenAtom::getNuclearAttractionEnergy() const
{
    return dft.getResult().nuclearAttractionEnergy;
}


double HydrogenAtom::getHartreeEnergy() const
{
    return dft.getResult().hartreeEnergy;
}


double HydrogenAtom::getExchangeEnergy() const
{
    return dft.getResult().exchangeEnergy;
}


// ============================================================
// DFT DATA
// ============================================================

const std::vector<double>&
HydrogenAtom::getRadialGrid() const
{
    return dft.getResult().r;
}


const std::vector<double>&
HydrogenAtom::getOrbital() const
{
    return dft.getResult().orbital;
}


const std::vector<double>&
HydrogenAtom::getDensity() const
{
    return dft.getResult().density;
}


const std::vector<double>&
HydrogenAtom::getEffectivePotential() const
{
    return dft.getResult().effectivePotential;
}


const DFTHydrogen::Result&
HydrogenAtom::getDFTResult() const
{
    return dft.getResult();
}