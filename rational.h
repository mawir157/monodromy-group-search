#pragma once

class Rational
{
  public:
    Rational(int m, unsigned int n);
    double value() const { return (double(m) / double(n)); }
    comp_d exp() const;
    bool operator<=(const Rational& r) const {return value() <= r.value();}
    bool operator<(const Rational& r) const {return value() < r.value();}
    std::string asString() const;

  public:
    int m;
    unsigned int n;
};

std::vector<Rational> GetFirstRationals(const size_t max,
                                        const bool skip_zero = true);

class Triple
{
  public:
    Triple(int m1, unsigned int n1, int m2, unsigned int n2,
           int m3, unsigned int n3);
    Triple(Rational r1, Rational r2, Rational r3);
    void reset(int m1, unsigned int n1, int m2, unsigned int n2,
                       int m3, unsigned int n3);
    double score() const {return r1.value() + r2.value() + r3.value();}
    Triple conj() const;
    Triple rot_by_pi() const;
    bool is_minimal() const;
    std::string asString() const;

  Rational r1;
  Rational r2;
  Rational r3;
};

bool MinimalSextuple(const Triple A, const Triple B);