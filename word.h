#pragma once
#include "includes.h"
#include "matrixfns.h"

class Word
{
  private:
    Word(const CompMat3& mat, const CompMat3& H);
    const CompMat3 mat;
    const CompMat3 H;
    const IsomClass cl;
    const int ord;
  public:
    CompMat3 matrix() const { return mat; }
    CompMat3 herm() const { return H; }
    int order() const { return ord; }
    IsomClass iso_class() const { return cl; }
    Word inv() const;

    bool operator==(const Word& rhs);
    Word operator*(const Word& rhs);
};