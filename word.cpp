#include "word.h"

Word::Word(const CompMat3& mat, const CompMat3& H) :
    mat(mat)
  , H(H)
  , cl(GetIsomClass(mat, H))
  , ord(Order(mat, H)) {}

Word Word::inv() const
{
  const CompMat3 m_inv = arma::inv(mat);
  const Word inv_mat(m_inv, H);
  return inv_mat;
}

bool Word::operator==(const Word& rhs)
{
  return areEqual(mat, rhs.matrix());
}

Word Word::operator*(const Word& rhs)
{
  const CompMat3 AB = mat * rhs.matrix();
  const Word prod_mat(AB, H);
  return prod_mat;
}