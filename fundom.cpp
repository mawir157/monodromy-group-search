#include "includes.h"
#include "fundom.h"
#include "matrixfns.h"

FunDom::FunDom(const Point& p, const CompMat3& H) :
    m_base(p)
  , m_H(H) {}


void FunDom::addPoints(const std::vector<Word>& wds)
{
  for (size_t i = 0; i <  wds.size(); ++i)
    this->addPoint(wds[i].get_matrix());
}

void FunDom::addPoint(const CompMat3& M)
{
  Point image = M * m_base;
  if (arma::norm(image - m_base) < TOL)
      return;
  for (size_t i = 0; i < m_images.size(); ++i)
  {
    // is already seem
    if (arma::norm(image - m_images[i]) < TOL)
      return;
  }
  m_images.push_back(image);
}

void FunDom::print() const
{
  for (size_t i = 0; i < m_images.size(); ++i)
  {
    if (isnan(comp_distance(m_images[i], m_base, m_H)))
    {
      std::cout << "**********************************************" << std::endl;
      std::cout << m_images[i] << std::endl;
      std::cout << m_base <<std::endl;
      const comp_d t = (herm(m_images[i], m_base, m_H) * herm(m_base, m_images[i], m_H)) /
                       (herm(m_images[i], m_images[i], m_H) * herm(m_base, m_base, m_H));
      std::cout << "<><> / <><> = " << std::real(t) << std::endl;
      std::cout << "sqrt(<><>/<><>)" << std::sqrt(std::real(t)) << std::endl;
    }
    std::cout << comp_distance(m_images[i], m_base, m_H) << std::endl;
  }
}

bool FunDom::stochastic_lattice(const unsigned int n,
                                const double width,
                                const bool verbose)
{
  // step 1 generate n random negative vectors
  std::vector<Point> neg_vectors;
  neg_vectors.push_back(m_base);

  double lower_bound = -width;
  double upper_bound = width;
  std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
  std::random_device rd;
  std::mt19937 gen(rd());

  while (neg_vectors.size() < n)
  {
    Point p(3);
    p << m_base[0] + comp_d(unif(gen), unif(gen)) << arma::endr
      << m_base[1] + comp_d(unif(gen), unif(gen)) << arma::endr
      << m_base[2] + comp_d(unif(gen), unif(gen)) << arma::endr;

    if (std::real(herm(p, p, m_H)) < 0)
      neg_vectors.push_back(p);
  }

  // now see how many of these points lie outside the dirichlet domain
  std::cout << "generated " << n << " negative vectors" << std::endl;

  std::vector<double> inside_distances;
  inside_distances.reserve(neg_vectors.size());
  for (unsigned int i = 0; i < neg_vectors.size(); ++i)
  {
    Point p = neg_vectors[i];
    bool ok = true;
    const double d = comp_distance(m_base, p, m_H);
    for (unsigned int j = 0; j < m_images.size(); ++j)
    {
      if (d > comp_distance(m_images[j], p, m_H))
      {
        ok = false;
        break;
      }
    }
    if (ok)
      inside_distances.push_back(d);
  }
  std::sort(inside_distances.begin(), inside_distances.end());

  if (verbose)
  {
    std::cout << " | " << inside_distances.size() << "~" << neg_vectors.size() << std::endl;
    std::cout << " | " << inside_distances.back() << std::endl;
  }
  return true;
}

std::vector<double> get_spectrum(const Point& p,
                                 const std::vector<Word>& words)
{
  std::vector<double> distances;
  distances.reserve(words.size());
  for (size_t i = 0; i <  words.size(); ++i)
  {
    const Word wd = words[i];
    const Point wd_p = wd.image(p);

    const double d = comp_distance(p, wd_p, wd.get_H_matrix());
    distances.push_back(d);
  }

  std::sort(distances.begin(), distances.end());

  std::vector<double> reduced_deistances;
  reduced_deistances.push_back(distances[0]);

  for (size_t i = 1; i < distances.size(); ++i)
  {
    if (std::abs(distances[i] - reduced_deistances.back()) > TOL)
      reduced_deistances.push_back(distances[i]);
  }

  return reduced_deistances;
}