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