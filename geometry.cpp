#include "geometry.h"

Polar normalize(const Point& p, size_t index)
{
  Polar norm(3);
  norm << p[0] / p[index] << arma::endr
       << p[1] / p[index] << arma::endr
       << p[2] / p[index] << arma::endr;

  return norm;
}

Polar get_polar(const Point& p1, const Point& p2, const CompMat3& H)
{
  arma::cx_drowvec lhs1 = p1.t() * H;
  arma::cx_drowvec lhs2 = p2.t() * H;

  comp_d a = lhs1[0];
  comp_d b = lhs1[1];
  comp_d c = lhs1[2];

  comp_d d = lhs2[0];
  comp_d e = lhs2[1];
  comp_d f = lhs2[2];

  comp_d z2 = ((c * d) - (a * f)) / ((a * e) - (b * d));
  comp_d z1 = -(c / a) - (b / a) * z2;

  Polar pol(3);
  pol << z1 << arma::endr
      << z2 << arma::endr
      << 1 << arma::endr;

  return pol;
}

Point intersect_C_lines(const Polar& p1, const Polar& p2, const CompMat3& H)
{

}

// intersection of B(*, 1), B(*, 2) and B(*, 3)
void triple_intersection(const Point& base, const Point& p1, const Point& p2,
                         const Point& p3, const CompMat3& H)
{
  // step 1 make sure we're generic i.e
  // no pair of 1, 2 and 3 are cospinal with 0
  Polar pol1 = get_polar(base, p1, H);
  Polar pol2 = get_polar(base, p2, H);
  Polar pol3 = get_polar(base, p3, H);

  // std::cout << "p0:" << std::endl
  //           << base << std::endl
  //           << "p1:" << std::endl
  //           << p1 << std::endl
  //           << "p2:" << std::endl
  //           << p2 << std::endl;

  CompMat3 P(3, 3);
  P << base[0] << p1[0] << p2[0] << arma::endr
    << base[1] << p1[1] << p2[1] << arma::endr
    << base[2] << p1[2] << p2[2] << arma::endr;

  comp_d u1 = herm(p3, p1, H) / herm(p3, base, H);
  comp_d u2 = herm(p3, p2, H) / herm(p3, base, H);

  // arma::cx_dvec e0(3);
  // e0 << 1 << arma::endr
  //    << 0 << arma::endr
  //    << 0 << arma::endr;

  // arma::cx_dvec e1(3);
  // e1 << 0 << arma::endr
  //    << 1 << arma::endr
  //    << 0 << arma::endr;

  // arma::cx_dvec e2(3);
  // e2 << 0 << arma::endr
  //    << 0 << arma::endr
  //    << 1 << arma::endr;

  // std::cout << "Pe0:" << std::endl;
  // std::cout << (P * e0) << std::endl;
  // std::cout << "Pe1:" << std::endl;
  // std::cout << (P * e1) << std::endl;
  // std::cout << "Pe2:" << std::endl;
  // std::cout << (P * e2) << std::endl;

  arma::cx_dvec U(3);
  U << 1 << arma::endr
    << u1 << arma::endr
    << u2 << arma::endr;

  CompMat3 A = arma::inv(P.t() * H);

  arma::cx_dvec test = A * U;

  std::cout << "U:" << std::endl;
  std::cout << (U) << std::endl;
  std::cout << "p3:" << std::endl;
  std::cout << normalize(p3) << std::endl;
  std::cout << "AU:" << std::endl;
  std::cout << normalize(test) << std::endl;

  std::cout << (P.t() * H * p3) << std::endl;
}
