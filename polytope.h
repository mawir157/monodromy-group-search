#pragma once
#include "includes.h"

class CWCell
{
  public:
    CWCell();
    explicit CWCell(const std::string name, const unsigned int id,
                    const std::vector<CWCell>& b, const unsigned int dim);
    const std::string name;
    const unsigned int id;
    const unsigned int dim;
    const std::vector<CWCell> boundary;
    int Euler() const;
};

void EulerRecursive(const std::vector<CWCell>& complex,
                    std::vector<std::set<unsigned int>>& cell_count);
