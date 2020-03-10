#pragma once
#include "includes.h"
#include "rational.h"

struct mat_sig
{
  mat_sig(unsigned int p, unsigned int n, unsigned int g);
  std::string asString() const;

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

CompMat3 TripleToMatrix(const Triple A, const bool norm=true);
CompMat3 SextupleToH(const Triple A, const Triple B, const bool force = true);
CompMat3 GetU(const Triple A);
double goldman(const comp_d& t);
bool isId(const CompMat3& matrix);
bool areEqual(const CompMat3& A, const CompMat3& B);
int Order(const CompMat3& m, const CompMat3& H,
          const int max_order = MAX_ORDER);
CompMat3 conj(const CompMat3 A, const CompMat3 P);
std::string GetIsomClassStr(const CompMat3& m, const CompMat3& H);
IsomClass GetIsomClass(const CompMat3& m, const CompMat3& H, const bool debug = false);
bool isPower(const CompMat3& A, const CompMat3& B,
             const unsigned int upto = 20);
int Braid(const CompMat3& A, const CompMat3& B,
          const unsigned int max_braid = MAX_BRAID);
mat_sig get_mat_sig(const CompMat3& matrix, const double tol=TOL);
bool is_non_finite_elliptic(const CompMat3& m, const CompMat3& H);
CompMat3 WordToMatrix(const unsigned int word,
                      const std::vector<CompMat3>& Mats);
bool CheckWords2(const unsigned int upto,
                 const std::vector<CompMat3>& Mats,
                 const CompMat3& H,
                 bool& SeenParabolic);
bool CheckWordsAlt(const std::vector<CompMat3>& Mats,
                   const CompMat3& H,
                   bool& SeenParabolic);
void getShortWords(const std::vector<CompMat3>& Generators,
                   std::vector<CompMat3>& Matrices,
                   const size_t upto = 1000);
std::string HumanWord(const unsigned int word,
                      const std::vector<std::string>& Names);
bool IsWordMinimal(const unsigned int word,
                   const size_t n);
void printShortWords(const unsigned int upto,
                     const std::vector<CompMat3>& Mats,
                     const std::vector<std::string>& Names,
                     const CompMat3& H,
                     const bool skip_loxo = true);
unsigned int StringToCode(const std::string word,
                          const std::vector<std::string>& Names);

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
void FindIsometryThatFixes(const Point& p,
                           const unsigned int upto,
                           const std::vector<CompMat3>& Mats,
                           const std::vector<std::string>& Names,
                           const CompMat3& H);