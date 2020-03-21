#pragma once
#include "geometry.h"
#include "matrixfns.h"
#include "word.h"

bool jorgensen(const std::vector<CompMat3>& words,
               const CompMat3& H,
               const bool verbose = false);
bool jorgensen_wrapper(const std::vector<CompMat3>& Mats,
                       const CompMat3& H,
                       const unsigned int upto,
                       const bool verbose = false);
bool jorgensen_all(const std::vector<CompMat3>& Mats,
                   const CompMat3& H,
                   const unsigned int upto);

bool jorgensen_new(const std::vector<Word>& Mats);