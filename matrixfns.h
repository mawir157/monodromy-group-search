#pragma once
#include "includes.h"
#include "rational.h"
#include "geometry.h"

struct mat_sig
{
  mat_sig(unsigned int p, unsigned int n, unsigned int g);
  std::string asString() const;
  bool match(unsigned int pos, unsigned int nul, unsigned int neg) const;

  unsigned int m_pos;
  unsigned int m_nul;
  unsigned int m_neg;
};

enum IsomClass { Identity,
                 Elliptic,
                 ReflectionPoint,
                 ReflectionLine,
                 ParabolicPure,
                 ParabolicScrew,
                 Loxodromic,
                 Failure };

enum Generator { ID, A, Ai, B, Bi, R1, R2, R3, E1, E2, E3 };
Generator inverse(const Generator gen);

CompMat3 TripleToMatrix(const Triple A, const bool norm=true);
CompMat3 SextupleToH(const Triple A, const Triple B, const bool force = true);
CompMat3 GetU(const Triple A);
double goldman(const comp_d& t);
bool isId(const CompMat3& matrix, const double tol=TOL);
bool areEqual(const CompMat3& A, const CompMat3& B, const double tol=TOL);
bool traceEqual(const comp_d& trace_a, const comp_d& trace_b, const double tol=TOL);
int base_order(const CompMat3& m, const int max_order = MAX_ORDER);
int Order(const CompMat3& m, const CompMat3& H,
          const int max_order = MAX_ORDER);
CompMat3 power(const CompMat3 A, const size_t p);
CompMat3 conj(const CompMat3 A, const CompMat3 P);
std::string GetIsomClassStr(const CompMat3& m, const CompMat3& H);
IsomClass GetIsomClass(const CompMat3& m, const CompMat3& H);
bool isPower(const CompMat3& A, const CompMat3& B,
             const unsigned int upto = 20);
int Braid(const CompMat3& A, const CompMat3& B,
          const unsigned int max_braid = MAX_BRAID);
mat_sig get_mat_sig(const CompMat3& matrix, const double tol=TOL);

comp_d herm(const Point& p1, const Point& p2, const CompMat3& H);
comp_d crossprod(const Point& z1, const Point& z2,
                 const Point& w1, const Point& w2, const CompMat3& H);
double Ma(const CompMat3& A);
double comp_distance(const Point& p1, const Point& p2, const CompMat3& H);
Point find_null(const CompMat3& H);
Point getEllipticFixedPoint(const CompMat3& M, const CompMat3& H);
Point getLineReflectionPolar(const CompMat3& M, const CompMat3& H);
void getLoxodromicFixed(const CompMat3& M, const CompMat3& H,
                        Point& p1, Point& p2);

void RemoveNearZeros(CompMat3& M, const double tol=TOL);
void RemoveNearIntegers(CompMat3& M, const double tol=TOL);
