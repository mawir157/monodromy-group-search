#pragma once

class FunDom
{
  public:
    FunDom(const Point& p, const CompMat3& H);
    const Point base;
    const CompMat3 H;
    std::vector<Point> images;
    void addPoint(CompMat3& M);
    void print() const;

    bool stochastic_lattice(unsigned int n, bool verbose = false);
};