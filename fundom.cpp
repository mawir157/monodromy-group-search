#include "includes.h"
#include "fundom.h"
#include "matrixfns.h"

FunDom::FunDom(const Point& p, const CompMat3& H) :
    base(p)
  , H(H) {}

void FunDom::addPoint(CompMat3& M)
{
  Point image = M * base;
  if (arma::norm(image - base) < TOL)
      return;
  for (size_t i = 0; i < images.size(); ++i)
  {
    // is already seem
    if (arma::norm(image - images[i]) < TOL)
      return;
  }
  images.push_back(image);
}

void FunDom::print() const
{
  for (size_t i = 0; i < images.size(); ++i)
  {
    if (isnan(comp_distance(images[i], base, H)))
    {
      std::cout << "**********************************************" << std::endl;
      std::cout << images[i] << std::endl;
      std::cout << base <<std::endl;
      const comp_d t = (herm(images[i], base, H) * herm(base, images[i], H)) / 
                   (herm(images[i], images[i], H) * herm(base, base, H));      
      std::cout << "<><> / <><> = " << std::real(t) << std::endl;
      std::cout << "sqrt(<><>/<><>)" << std::sqrt(std::real(t)) << std::endl;
    }
    std::cout << comp_distance(images[i], base, H) << std::endl;
  }
}

bool FunDom::stochastic_lattice(unsigned int n, bool verbose)
{
  // step 1 generate n random negative vectors
  std::vector<Point> null_vectors;

  double lower_bound = -10;
  double upper_bound = 10; 
  std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
  std::random_device rd;
  std::mt19937 gen(rd());

  while (null_vectors.size() < n)
  {
    Point p(3);
    p << base[0] + comp_d(unif(gen), unif(gen)) << arma::endr
      << base[1] + comp_d(unif(gen), unif(gen)) << arma::endr
      << base[2] + comp_d(unif(gen), unif(gen)) << arma::endr;

    if (std::real(herm(p, p, H)) < 0)
      null_vectors.push_back(p);
  }

  // now see how many of these points lie outside the dirichlet domain
  unsigned int c = 0;
  for (unsigned int i = 0; i < null_vectors.size(); ++i)
  {
    Point p = null_vectors[i];
    for (unsigned int j = 0; j < images.size(); ++j)
    {
      if (comp_distance(base, p, H) > comp_distance(images[j], p, H))
      {
        ++c;
        break;
      }
    }
  }
  if (verbose)
    std::cout << " | " << c  << "~" << null_vectors.size() << std::endl;
  return true;
}