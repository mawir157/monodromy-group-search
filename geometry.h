#pragma once

#include "includes.h"
#include "matrixfns.h"

Polar normalize(const Point& p1, size_t index = 0);
Polar get_polar(const Point& p1, const Point& p2, const CompMat3& H);
Point intersect_C_lines(const Polar& p1, const Polar& p2, const CompMat3& H);
void triple_intersection(const Point& base, const Point& p1, const Point& p2,
                         const Point& p3, const CompMat3& H);
