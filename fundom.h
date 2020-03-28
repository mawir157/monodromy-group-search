#pragma once

#include "word.h"

class FunDom
{
  public:
    FunDom(const Point& p, const CompMat3& H);

    const Point        m_base;
    const CompMat3     m_H;
    std::vector<Point> m_images;

    void addPoints(const std::vector<Word>& wds);
    void print() const;
    bool stochastic_lattice(const unsigned int n,
														const double width = 100.0,
														const bool verbose = false);

  private:
    void addPoint(const CompMat3& M);
};

std::vector<double> get_spectrum(const Point& p,
	                               const std::vector<Word>& words);
