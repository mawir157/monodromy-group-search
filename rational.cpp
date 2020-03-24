#include "includes.h"
#include "utils.h"
#include "rational.h"

Rational::Rational(int m, unsigned int n) :
    m(m)
  , n(n) {}

std::string Rational::asString(const bool formal) const
{
  if ((m == 0) && (formal))
    return "0";

  std::string s = std::to_string(m);
  if (formal)
    s.append(",");
  else
    s.append("/");

  s.append(std::to_string(n));

  return s;
}

comp_d Rational::exp() const
{
   const comp_d theta(0.0, 2.0 * PI * value());
   return std::exp(theta);
}

std::vector<Rational> GetFirstRationals(const size_t max,
                                        const bool skip_zero)
{
  std::vector<Rational> rs;
  unsigned int init = (skip_zero) ? 1 : 0;
  for(unsigned int im = init; im <= max; ++im)
  {
    for (unsigned int in = 1; in <= max; ++in)
    {
      if (im >= in)
        continue;

      if (gcd(im, in) == 1)
        rs.emplace_back(im, in);

      if (im == 0)
        break;
    }
  }

  return rs;
}

Triple::Triple(int m1, unsigned int n1, int m2, unsigned int n2,
               int m3, unsigned int n3) :
    r1(Rational(m1, n1))
  , r2(Rational(m2, n2))
  , r3(Rational(m3, n3)) {}

Triple::Triple(Rational r1, Rational r2, Rational r3) :
    r1(r1)
  , r2(r2)
  , r3(r3) {}

void Triple::reset(int m1, unsigned int n1, int m2, unsigned int n2,
                   int m3, unsigned int n3)
{
  r1 = Rational(m1, n1);
  r2 = Rational(m2, n2);
  r3 = Rational(m3, n3);
}

Triple Triple::conj() const
{
  return Triple(r1.m == 0 ? 0 : r1.n - r1.m, r1.n,
                r2.m == 0 ? 0 : r2.n - r2.m, r2.n,
                r3.m == 0 ? 0 : r3.n - r3.m, r3.n);
}

Triple Triple::rot_by_pi() const
{
  // add half to each rational
  int new_m1 = 2 * r1.m + r1.n;
  unsigned int new_n1 = 2 * r1.n;
  if (2 * r1.m + r1.n >= 2 * r1.n)
    new_m1 -= 2 * r1.n;
  if ((new_m1 % 2 == 0))
  {
    new_m1 /= 2;
    new_n1 /= 2;
  }

  int new_m2 = 2 * r2.m + r2.n;
  unsigned int new_n2 = 2 * r2.n;
  if (2 * r2.m + r2.n >= 2 * r2.n)
    new_m2 -= 2 * r2.n;
  if ((new_m2 % 2 == 0))
  {
    new_m2 /= 2;
    new_n2 /= 2;
  }

  int new_m3 = 2 * r3.m + r3.n;
  unsigned int new_n3 = 2 * r3.n;
  if (2 * r3.m + r3.n >= 2 * r3.n)
    new_m3 -= 2 * r3.n;
  if ((new_m3 % 2 == 0))
  {
    new_m3 /= 2;
    new_n3 /= 2;
  }

  return Triple(new_m1, new_n1,
                new_m2, new_n2,
                new_m3, new_n3);
}

bool Triple::is_minimal() const
{
  Triple c = conj();
  Triple r = rot_by_pi();
  Triple rc = r.conj();

  return ((score() <= c.score()) &&
          (score() <= r.score()) &&
          (score() <= rc.score()));
}

std::string Triple::asString() const
{
  std::string s = "(";
  s.append(r1.asString());
  s.append(", ");
  s.append(r2.asString());
  s.append(", ");
  s.append(r3.asString());
  s.append(")");

  return s;
}

std::string Triple::formal() const
{
  std::string s = "";
  s.append(r1.asString(true));
  s.append(",");
  s.append(r2.asString(true));
  s.append(",");
  s.append(r3.asString(true));
  s.append(",");

  return s;
}

bool MinimalSextuple(const Triple A, const Triple B)
{
  const double a = A.score();
  const double b = B.score();

  const double ar = A.rot_by_pi().score();
  const double br = B.rot_by_pi().score();

  const double ac = A.conj().score();
  const double bc = B.conj().score();

  const double arc = A.rot_by_pi().conj().score();
  const double brc = B.rot_by_pi().conj().score();

  return (((a + b) <= (ar + br)) &&
          ((a + b) <= (ac + bc)) &&
          ((a + b) <= (arc + brc)));
}
