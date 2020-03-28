#pragma once

#include "word.h"

class FunDom
{
  public:
    FunDom(const Point& p, const CompMat3& H);

    void addPoints(const std::vector<Word>& wds);
    void print() const;
    unsigned int stochastic_lattice(const unsigned int n,
													      	  const double width = 100.0,
														        const bool verbose = false) const ;
    std::vector<double> repel_lattice(const unsigned int n,
														          const unsigned int steps = 1000,
																			const bool verbose = false);

  private:
    const Point        m_base;
    const CompMat3     m_H;
    std::vector<Point> m_images;

    void addPoint(const CompMat3& M);
    void move_point(Point& p, const bool verbose = false);
    bool is_in_Dirichlet(const Point& p) const;
};

std::vector<double> get_spectrum(const Point& p,
	                               const std::vector<Word>& words);
