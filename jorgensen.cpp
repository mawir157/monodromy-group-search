#include "jorgensen.h"

bool jorgensen(const std::vector<Word>& Words)
{
  int count = 0;
  for (unsigned int i = 0; i < Words.size(); ++i)
  {
    const Word A = Words[i];

    for (unsigned int j = 0; j < Words.size(); ++j)
    {
      if (i == j)
        continue;

      const Word B = Words[j];
      const CompMat3 Am = A.get_matrix();
      const CompMat3 Bm = B.get_matrix();
      const CompMat3 H  = A.get_H_matrix();
      // If A is loxomodromic
      if (A.get_isom_class() == IsomClass::Loxodromic)
      {
        Point mu;
        Point nu;
        getLoxodromicFixed(Am, H, mu, nu);
        const double M_A = Ma(Am);

        if (M_A < 1)
        {
          ++count;
          // test 1
          comp_d cross = crossprod(Bm * mu, nu, mu, Bm * nu, H);
          const double test_1 = M_A * (std::sqrt(std::abs(cross)) + 1);
          if (test_1 < (1 + TOL))
          {
            std::cout << "Faliure 1" << std::endl;
            std::cout << A.as_string() << " ~ ";
            std::cout << B.as_string() << " = ";
            std::cout << test_1 - 1 << std::endl;
            std::cout << A.get_matrix() << std::endl;
            std::cout << arma::trace(A.get_matrix()) << std::endl;
            std::cout << goldman(arma::trace(A.get_matrix())) << std::endl;
            return false;
          }

          // test 2
          cross = crossprod(Bm * mu, mu, nu, Bm * nu, H);
          const double test_2 = M_A * (std::sqrt(std::abs(cross)) + 1);
          if (test_2 < (1 + TOL))
          {
            std::cout << "Faliure 2" << std::endl;
            return false;
          }

          // test 3
          // Not implemented

          // test 4
          cross = crossprod(mu, nu, Bm * mu, Bm * nu, H);
          const double test_4 = M_A + std::sqrt(std::abs(cross));
          if (test_4 < (1 + TOL))
          {
            std::cout << "Faliure 4" << std::endl;
            return false;
          }
        }
        //std::cout << "(" << test_1  << "," << test_2 << "," << test_4 << ")" << std::endl;
      }

      else if (A.get_isom_class() == IsomClass::ReflectionLine)
      {

      }

      else if (A.get_isom_class() == IsomClass::Elliptic)
      {
        if (A.get_order() < 1) // this should never be hit
          return false; // non-finite-order-elliptic => non-discrete

        if (A.get_order() <= 7) // can't do jorgensen if order is less than 7
          continue;

        ++count;
        Point z = getEllipticFixedPoint(Am, H);
        Point new_z = normalize(Bm * z);

        // If B fixes z can't uses jorgensen
        if (arma::approx_equal(z, new_z, "absdiff", TOL))
          continue;

        double jorg_value = std::cosh(comp_distance(new_z, z, H) / 2) *
                            std::sin(PI / 2.0);

        if (jorg_value <  0.5)
        {
          std::cout << "Faliure Elliptic" << std::endl;
          return false;
        }
      }
    }
  }
  std::cout << "[[" << count << "]]" << std::endl;
  return true;
}
