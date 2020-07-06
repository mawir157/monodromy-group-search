#pragma once

#include "includes.h"
#include "word.h"

class StochFunDom
{
  public:
    StochFunDom(const Point& p, const CompMat3& H);

    void addPoints(const std::vector<Word>& wds);
    void print() const;
    unsigned int stochastic_lattice(const unsigned int n,
						            const double width = 100.0,
								    const bool verbose = false) const ;
    std::vector<double> repel_lattice(const unsigned int n,
								      const unsigned int steps = 1000,
								      const bool verbose = false) const;
    void fixed_points() const;

  private:
    const Point        m_base;
    const CompMat3     m_H;
    std::vector<Point> m_images;

    void addPoint(const CompMat3& M);
    void move_point(Point& p, const bool verbose = false) const;
    bool is_in_Dirichlet(const Point& p) const;

};

std::vector<double> get_spectrum(const Point& p,
	                               const std::vector<Word>& words);

std::vector<Point> get_orbit(const Point& p,
                             const std::vector<Word>& words);

class FunDom
{
  public:
    FunDom(const std::vector<Word>& words, const CompMat3& H);
    std::vector<Word> get_vertices(const bool verbose = false) const;

  private:
    const std::vector<Word> m_words;
    const CompMat3              m_H;
};