#include "jorgensen.h"

bool jorgensen(const std::vector<CompMat3>& words, const CompMat3& H,
               const bool verbose)
{
  bool any = false;
  for (unsigned int i = 0; i < words.size(); ++i)
  {
    const CompMat3 A = words[i];
    const comp_d tr = arma::trace(A);

    // check A has real trace
    if (std::abs(std::imag(tr)) > TOL)
    {
      //std::cout << "t" << std::endl;
      continue;
    }

    const IsomClass cl = GetIsomClass(A, H);
    // check A is regular ellitpic
    if (cl != IsomClass::Elliptic)
    {
      //std::cout << "e" << std::endl;
      continue;
    }

    // check A has finite order
    const int n = Order(A, H);
    if (n < 1)
      return false;

    if (n <= 7)
    {
      //std::cout << "7" << std::endl;
      continue;
    }

    any |= true;

    Point z = getEllipticFixedPoint(A, H);

    // now find a word that does not fix z
    for (unsigned int j = 0; j < words.size(); ++j)
    {
      if (i == j)
        continue;

      const CompMat3 B = words[j];
      Point new_z = normalize(B * z);

      // B fixes z
      if (arma::approx_equal(z, new_z, "absdiff", TOL))
        continue;

      double jorg_value = std::cosh(comp_distance(new_z, z, H) / 2) *
                          std::sin(PI / 2.0);

      if (jorg_value <  0.5)
        return false;
    }
  }
  if (!any && verbose)
    std::cout << "Got to the end without ever applying J-test :(" << std::endl;
  return true;
}

bool jorgensen_wrapper(const std::vector<CompMat3>& Mats,
                       const CompMat3& H,
                       const unsigned int upto,
                       const bool verbose)
{
  std::vector<CompMat3> matrices;
  for (unsigned int i = 0; i <= upto; ++i)
  {
    bool minimal = IsWordMinimal(i, Mats.size());
    if (!IsWordMinimal(i, Mats.size()))
       continue;

    const CompMat3 A = WordToMatrix(i, Mats);
    bool seen = false;
    int j = 0;
    for (j = 0; j < matrices.size(); ++j)
    {
      const CompMat3 B = matrices[j];
      seen = areEqual(A, B);
      if (seen)
        break;
    }

    if (!seen)
      matrices.push_back(A);
  }

  bool j_test = jorgensen(matrices, H, verbose);
  if (verbose)
  {
    std::cout << matrices.size() << std::endl;
    std::cout << (j_test ? "PASSED" : "FAILED") << std::endl;
  }
  return j_test;
}

bool jorgensen_all(const std::vector<CompMat3>& Mats,
                   const CompMat3& H,
                   const unsigned int upto)
{
  bool any = false;
  std::vector<CompMat3> matrices;
  for (unsigned int i = 0; i <= upto; ++i)
  {
    bool minimal = IsWordMinimal(i, Mats.size());
    if (!IsWordMinimal(i, Mats.size()))
       continue;

    const CompMat3 A = WordToMatrix(i, Mats);
    bool seen = false;
    int j = 0;
    for (j = 0; j < matrices.size(); ++j)
    {
      const CompMat3 B = matrices[j];
      seen = areEqual(A, B);
      if (seen)
        break;
    }

    if (!seen)
      matrices.push_back(A);
  }

  for (unsigned int i = 0; i < matrices.size(); ++i)
  {
    const CompMat3 A = matrices[i];
    const IsomClass cl_A = GetIsomClass(A, H);

    for (unsigned int j = 0; j < matrices.size(); ++j)
    {
      const CompMat3 B = matrices[j];

      if (cl_A == IsomClass::Loxodromic)
      {
        Point mu;
        Point nu;
        getLoxodromicFixed(A, H, mu, nu);
        const double M_A = Ma(A);

        if (M_A < 1)
        {
          // test 1
          comp_d cross = crossprod(B * mu, nu, mu, B * nu, H);
          const double test_1 = M_A * (std::sqrt(std::abs(cross)) + 1);
          if (test_1 < (1 + TOL))
            return false;

          // test 2
          cross = crossprod(B * mu, mu, nu, B * nu, H);
          const double test_2 = M_A * (std::sqrt(std::abs(cross)) + 1);
          if (test_2 < (1 + TOL))
            return false;

          // test 3
          // Not implemented

          // test 4
          cross = crossprod(mu, nu, B * mu, B * nu, H);
          const double test_4 = M_A + std::sqrt(std::abs(cross));
          if (test_4 < (1 + TOL))
            return false;
        }
        //std::cout << "(" << test_1  << "," << test_2 << "," << test_4 << ")" << std::endl;
      }

      if (cl_A == IsomClass::ReflectionLine)
      {

      }
    }
  }
  return true;
}

bool jorgensen_new(const std::vector<Word>& Words)
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
