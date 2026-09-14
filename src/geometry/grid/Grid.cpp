#include "Grid.h"

#include "GeometryConstants.h"


Grid::Grid()
    : size(GeometryConstants::GRID_SIZE),
      divisions(GeometryConstants::GRID_DIVISIONS),
      spacing(GeometryConstants::GRID_SPACING),
      min(-size * 0.5f),
      max(size * 0.5f)
{
}


float Grid::getSize() const
{
    return size;
}


int Grid::getDivisions() const
{
    return divisions;
}


float Grid::getSpacing() const
{
    return spacing;
}


const glm::vec3& Grid::getMin() const
{
    return min;
}


const glm::vec3& Grid::getMax() const
{
    return max;
}