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

unsigned int FunDom::stochastic_lattice(const unsigned int n,
                                        const double width,
                                        const bool verbose) const
{
  // step 1 generate n random negative vectors
  std::vector<Point> neg_vectors;
  neg_vectors.push_back(m_base);

  std::uniform_real_distribution<double> unif(-width, width);
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
    const double d = comp_distance(m_base, p, m_H);
    if (this->is_in_Dirichlet(p))
      inside_distances.push_back(d);
  }
  std::sort(inside_distances.begin(), inside_distances.end());

  if (verbose)
  {
    std::cout << " | " << inside_distances.size() << "~" << neg_vectors.size() << std::endl;
    std::cout << " | " << inside_distances.back() << std::endl;
  }
  return inside_distances.size();
}

std::vector<double> FunDom::repel_lattice(const unsigned int n,
                                          const unsigned int steps,
                                          const bool verbose)
{
  std::vector<double> distances;

  for (size_t i = 0; i < n; ++i) 
  {
    Point p = find_neg(m_H);
    while (!is_in_Dirichlet(p))
      p = find_neg(m_H);
    
    std::cout << i << std::endl;
    for (size_t j = 0; j < steps; ++j)
    {
      move_point(p, verbose);
    }

    const double d = comp_distance(p, m_base, m_H);
    distances.push_back(d);    
  }

  std::sort(distances.begin(), distances.end());
  if (verbose)
    std::cout << " ------> " << distances.back() << std::endl;

  return distances;
}

void FunDom::move_point(Point& p, const bool verbose)
{
  std::uniform_real_distribution<double> unif(-0.1, 0.1);
  std::random_device rd;
  std::mt19937 gen(rd());

  const double dist = comp_distance(p, m_base, m_H);
if (verbose) {std::cout << "]--> " << dist << std::endl;}

  while (true)
  {
    Point new_p(3);
    new_p << p[0] + comp_d(unif(gen), unif(gen)) << arma::endr
          << p[1] + comp_d(unif(gen), unif(gen)) << arma::endr
          << p[2] + comp_d(unif(gen), unif(gen)) << arma::endr;
    if (std::real(herm(new_p, new_p, m_H)) >= 0)
      continue;

    double new_dist = comp_distance(new_p, m_base, m_H);
    if (new_dist <= dist)
      continue;

    if(!is_in_Dirichlet(new_p))
      continue;

if (verbose) {std::cout << "|---------------->" << new_dist << std::endl;}
    p = new_p;
    break;
  }
}

bool FunDom::is_in_Dirichlet(const Point& p) const
{
  bool ok = true;
  const double d = comp_distance(m_base, p, m_H);
  for (unsigned int j = 0; j < m_images.size(); ++j)
    if (d > comp_distance(m_images[j], p, m_H))
      return false;

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