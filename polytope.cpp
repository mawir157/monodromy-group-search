#include "polytope.h"

CWCell::CWCell() :
    name("NULL")
  , id(1729)
  , dim(691) {};

CWCell::CWCell(const std::string name, const unsigned int id,
               const std::vector<CWCell>& b,const unsigned int dim) :
    name(name)
  , id(id)
  , dim(dim) 
  , boundary(b) {};

int CWCell::Euler() const
{
  std::vector<std::set<unsigned int>> cell_count(dim + 1);
  cell_count[dim].insert(id);

  for (size_t i = 0; i < boundary.size(); ++i)
  {
    const CWCell complex = boundary[i];
    cell_count[complex.dim].insert(complex.id);
    if (complex.dim != 0) 
      EulerRecursive(complex.boundary, cell_count);
  }

  int euler = 0;
  for (size_t i = 0; i < cell_count.size(); ++i)
    if (i % 2 ==0)
      euler += cell_count[i].size();
    else
      euler -= cell_count[i].size();

  return euler;
}

void EulerRecursive(const std::vector<CWCell>& boundary,
                    std::vector<std::set<unsigned int>>& cell_count)
{
  for (size_t i = 0; i < boundary.size(); ++i)
  {
    const CWCell complex = boundary[i];
    cell_count[complex.dim].insert(complex.id);
    if (complex.dim != 0) 
      EulerRecursive(complex.boundary, cell_count);
  }
}