#include "includes.h"
#include "rational.h"
#include "matrixfns.h"
#include "geometry.h"

CompMat3 TripleToMatrix(const Triple T,
                        const bool norm)
{
  const comp_d a1 = T.r1.exp();
  const comp_d a2 = T.r2.exp();
  const comp_d a3 = T.r3.exp();

  CompMat3 mat(3, 3);
  mat << a1 + a2 + a3 << -(a1 * a2 + a2 * a3 + a3 * a1) << a1 * a2 * a3 << arma::endr
      << 1.0          << 0.0                            << 0.0          << arma::endr
      << 0.0          << 1.0                            << 0.0          << arma::endr;

  if (norm)
  {
    comp_d lambda = arma::det(mat);
    lambda = std::pow(lambda, -1.0/3.0);
    mat = lambda * mat;
  }

  return mat;
}

CompMat3 GetU(const Triple A)
{
  const comp_d a1 = A.r1.exp();
  const comp_d a2 = A.r2.exp();
  const comp_d a3 = A.r3.exp();

  CompMat3 ubar(3, 3);
  ubar << 1.0 << -(a2 + a3) << a2 * a3 << arma::endr
       << 1.0 << -(a3 + a1) << a3 * a1 << arma::endr
       << 1.0 << -(a1 + a2) << a1 * a2 << arma::endr;

  return arma::inv(ubar);
}

CompMat3 SextupleToH(const Triple A, const Triple B, const bool force)
{
  const comp_d d1 = 2 * std::sin(PI * (B.r1.value() - A.r1.value())) *
    (std::sin(PI * (B.r2.value() - A.r1.value())) /
                                    std::sin(PI * (A.r2.value() - A.r1.value()))) *
    (std::sin(PI * (B.r3.value() - A.r1.value())) /
                                    std::sin(PI * (A.r3.value() - A.r1.value())));

  const comp_d d2 = 2 * std::sin(PI * (B.r2.value() - A.r2.value())) *
    (std::sin(PI * (B.r1.value() - A.r2.value())) /
                                    std::sin(PI * (A.r1.value() - A.r2.value()))) *
    (std::sin(PI * (B.r3.value() - A.r2.value())) /
                                    std::sin(PI * (A.r3.value() - A.r2.value())));

  const comp_d d3 = 2 * std::sin(PI * (B.r3.value() - A.r3.value())) *
    (std::sin(PI * (B.r2.value() - A.r3.value())) /
                                    std::sin(PI * (A.r2.value() - A.r3.value()))) *
    (std::sin(PI * (B.r1.value() - A.r3.value())) /
                                    std::sin(PI * (A.r1.value() - A.r3.value())));

  CompMat3 mat(3, 3);
  mat << d1  << 0.0 << 0.0 << arma::endr
      << 0.0 << d2  << 0.0 << arma::endr
      << 0.0 << 0.0 << d3  << arma::endr;

  CompMat3 u = GetU(A);

  mat = arma::inv(u.t()) * mat * arma::inv(u);

  // force the matrix to have sig (3,0) or (2,1) if possible
  if (force)
  {
    mat_sig sig = get_mat_sig(mat);
    if ((sig.m_pos == 0 && sig.m_neg == 3) || (sig.m_pos == 1 && sig.m_neg == 2))
      mat *= -1.0;
  }
  return mat;
}

double goldman(const comp_d& t)
{
  double a = std::abs(t);

  return (a * a * a * a) - (8 * std::real(t * t * t)) + (18 * a * a) - 27;
}

bool isId(const CompMat3& matrix)
{
  CompMat3 id = arma::eye<CompMat3>(3,3);
  if (arma::approx_equal(matrix, id, "absdiff", TOL))
    return true;

  if (arma::approx_equal(matrix, omega * id, "absdiff", TOL))
    return true;

  if (arma::approx_equal(matrix, omega * omega * id, "absdiff", TOL))
    return true;

  return false;
}

bool areEqual(const CompMat3& A, const CompMat3& B)
{
  if (arma::approx_equal(A, B, "absdiff", TOL))
    return true;

  if (arma::approx_equal(A, omega * B, "absdiff", TOL))
    return true;

  if (arma::approx_equal(A, omega * omega * B, "absdiff", TOL))
    return true;

  return false;  
}

int Order(const CompMat3& m, const CompMat3& H, const int max_order)
{
  if (isId(m))
    return 1;

  IsomClass cl = GetIsomClass(m, H);
  if ((cl != IsomClass::Elliptic) &&
      (cl != IsomClass::ReflectionPoint) &&
      (cl != IsomClass::ReflectionLine))
    return -1;

  CompMat3 m_copy = m;

  for (int i = 2; i < max_order; ++i)
  {
    m_copy = m_copy * m;
    if (isId(m_copy))
      return i;
  }

  return -2;
}

CompMat3 conj(const CompMat3 A, const CompMat3 P)
{
  return (P * A * arma::inv(P));
}

std::string GetIsomClassStr(const CompMat3& m, const CompMat3& H)
{
  const IsomClass isom = GetIsomClass(m, H);

  switch(isom)
  {
    case IsomClass::Elliptic:
    {
      std::string reg_ell = "Regular Elliptic";
      reg_ell.append(" <");
      reg_ell.append(std::to_string(Order(m, H)));
      reg_ell.append(">");
      return reg_ell;
      break;
    }
    case IsomClass::ReflectionPoint:
    {
      std::string reg_refp = "Reflection point";
      reg_refp.append(" <");
      reg_refp.append(std::to_string(Order(m, H)));
      reg_refp.append(">");
      return reg_refp;
      break;
    }
    case IsomClass::ReflectionLine:
    {
      std::string reg_refl = "Reflection line";
      reg_refl.append(" <");
      reg_refl.append(std::to_string(Order(m, H)));
      reg_refl.append(">");
      return reg_refl;
      break;
    }
    case IsomClass::ParabolicPure:
      return "Pure Parabolic";
      break;
    case IsomClass::ParabolicScrew:
      return "Screw Parabolic";
      break; 
    case IsomClass::Loxodromic:
      return "Loxodromic";
      break;
  case IsomClass::Identity:
      return "Identity";
      break;
    case IsomClass::Failure:
    default:
      return "Failure";
      break;
  }

  return "ERROR!";
}

IsomClass GetIsomClass(const CompMat3& m, const CompMat3& H, const bool debug)
{
// if (debug)
//   std::cout << "m = " << std::endl << m << std::endl;

// if (debug)
//   std::cout << "order = " << Order(m, H) << std::endl;

  if (isId(m))
    return IsomClass::Identity;

  comp_d tr = arma::trace(m);
// if (debug)
//   std::cout << "trace = " << tr << std::endl;
  if ((std::abs(tr - 3.0) < TOL) ||
      (std::abs(tr - 3.0 * omega) < TOL) ||
      (std::abs(tr - 3.0 * omega * omega) < TOL))
    return IsomClass::ParabolicPure;

  double gman = goldman(tr);
// if (debug)
//   std::cout << "goldman = " << gman << std::endl;
  if (gman > LOWER_TOL)
    return IsomClass::Loxodromic;

  if (gman < -LOWER_TOL)
    return IsomClass::Elliptic;

  // gman is "0" we're either screw parabolic or a complex reflection
  // the trace lies on the smooth part of the deltoid
  // this means the e-vals are A A B
  if (std::abs(gman) < LOWER_TOL)
  {
    arma::cx_vec eigval;
    arma::cx_mat eigvec;

    arma::eig_gen(eigval, eigvec, m);

// if (debug)
// {
//   std::cout << "e-vals" << std::endl << eigval << std::endl;
//   std::cout << std::abs(eigval[0] - eigval[1]) << std::endl;
//   std::cout << std::abs(eigval[1] - eigval[2]) << std::endl;
//   std::cout << std::abs(eigval[2] - eigval[0]) << std::endl;
// }

// if (debug)
//   std::cout << "e-vecs" << std::endl << eigvec << std::endl;

    int i_1 = -1;
    int i_2 = -1;
    int j = - 1;
    arma::cx_dvec v_pos(3);

    // we have to be a lot more lax about tolerances here
    // because calculating eigen values is difficult
    if (std::abs(eigval[0] - eigval[1]) < 1e-4)
    {
      i_1 = 0; i_2 = 1; j = 2;
    }
    else if (std::abs(eigval[1] - eigval[2]) < 1e-4)
    {
      i_1 = 1; i_2 = 2; j = 0;
    }
    else if (std::abs(eigval[2] - eigval[0]) < 1e-4)
    {
      i_1 = 2; i_2 = 0; j = 1;
    }
    else
    {
      std::cout << std::endl << "GetIsomClass has failed!" << std::endl;
      std::cout << "trace = " << std::endl << tr << std::endl;
      std::cout << "goldman = " << std::endl << gman << std::endl;
      std::cout << "e-vals = " << std::endl << eigval << std::endl;
      std::cout << "cross = " << std::endl << eigvec.t() * H * eigvec << std::endl;

      return IsomClass::Failure;
    }      

    CompMat3 mHm = eigvec.t() * H * eigvec;

// if (debug)
//   std::cout << "pairwise" << std::endl << mHm << std::endl;

    const double r1 = std::real(mHm(i_1, i_1));
    const double r2 = std::real(mHm(i_2, i_2));
    const double s = std::real(mHm(j, j));

// if (debug)
//   std::cout << "r1 & r2 = " << r1 << " & " << r2 << std::endl;

    if ((std::abs(r1) < TOL) || (std::abs(r2) < TOL))
    {
      if (s > TOL)
        return IsomClass::ParabolicScrew;
      else
      {
        std::cout << "Not a parabolic!" << std::endl;
        return IsomClass::Failure;
      }

    }

    // if the e-vectors have opposite signs then we're definitely
    // boundary ellitpic
    if (r1 * r2 < 0)
      return IsomClass::ReflectionLine;

    // at this point the e-vectors have the same sign
    // and now need to check the cross terms
    const comp_d rcross = mHm(i_1, i_2);

    // if both e-vecs are +ve and the cross term is zero
    // then the e-space is +ve and we're a relfection in a point
    // if the cross term is non-zero then we're a relfection in a line
    if ((r1 > TOL) && (r2 > TOL))
    {
      if (std::abs(rcross) < TOL)
        return IsomClass::ReflectionPoint;
      else
        return IsomClass::ReflectionLine;
    }

    // At this point both e-vectors are -ve
    // If the cross term is zero then the e-space is purely -ve
    // which feels like an error, if the cross term is non-zero
    // then we're a relfection in a line
    if (std::abs(rcross) > TOL)
      return IsomClass::ReflectionLine;

    std::cout << "We are about explode!" << std::endl;
    std::cout << "evals = " << std::endl << eigval << std::endl;
    std::cout << "cross matrix = " << std::endl << mHm << std::endl;

    return IsomClass::Failure;
  }

  arma::cx_vec eigval;
  arma::cx_mat eigvec;

  arma::eig_gen(eigval, eigvec, m);

  std::cout << std::endl << "GetIsomClass has completely failed!" << std::endl;
  std::cout << tr << std::endl;
  std::cout << gman << std::endl;
  std::cout << m << std::endl;
  std::cout << eigval << std::endl;
  std::cout << eigvec.t() * H * eigvec << std::endl;
  std::cout << std::abs(eigval[0] - eigval[1]) << std::endl;
  std::cout << std::abs(eigval[1] - eigval[2]) << std::endl;
  std::cout << std::abs(eigval[2] - eigval[0]) << std::endl;
  return IsomClass::Failure;
}

// is A = B^k or for some k positive or negative
bool isPower(const CompMat3& A, const CompMat3& B, const unsigned int upto)
{
  CompMat3 test = B;
  CompMat3 test_inv = arma::inv(B);

  for (unsigned int i = 1; i < upto; ++i)
  {
    if (areEqual(A, test))
      return true;

    if (areEqual(A, test_inv))
      return true;

    test *= B;
    test_inv *= arma::inv(B);
  }
  return false;
}

int Braid(const CompMat3& A, const CompMat3& B, const unsigned int max_braid)
{
  CompMat3 lhs = A;
  CompMat3 rhs = B;
  int braid = 1;
  bool parity = true;

  while (braid < max_braid)
  {
    if (areEqual(lhs, rhs))
      return braid;

    if (parity)
    {
      lhs *= B;
      rhs *= A;
      parity = false;
    }
    else
    {
      lhs *= A;
      rhs *= B;
      parity = true;
    }
    ++braid;
  }
  return -1;
}

mat_sig::mat_sig(unsigned int p, unsigned int n, unsigned int g) :
    m_pos(p)
  , m_nul(n)
  , m_neg(g) {}

std::string mat_sig::asString() const
{
  std::string s = "(";
  s.append(std::to_string(m_pos));
  s.append(", ");
  s.append(std::to_string(m_neg));
  s.append(")");

  return s;
}

mat_sig get_mat_sig(const CompMat3& matrix, const double tol)
{
  arma::cx_vec eigvals = arma::eig_gen( matrix );
  unsigned int pos_evals = 0;
  unsigned int nul_evals = 0;
  unsigned int neg_evals = 0;
  for (size_t i = 0; i < eigvals.size(); ++i)
  {
    if (std::abs(std::real(eigvals[i])) > TOL)
    {
      if (std::real(eigvals[i]) > 0.0)
        ++pos_evals;
      else
        ++neg_evals;
    }
    else
      ++nul_evals;
  }
  mat_sig ret(pos_evals, nul_evals, neg_evals);
  return ret;
}

bool is_non_finite_elliptic(const CompMat3& m, const CompMat3& H)
{
  const IsomClass cl = GetIsomClass(m, H);
  if ((cl== IsomClass::Elliptic) ||
      (cl == IsomClass::ReflectionPoint) ||
      (cl == IsomClass::ReflectionLine))
    return (Order(m, H) < 0);

  return false;
}

comp_d herm(const Point& p1, const Point& p2, const CompMat3& H)
{
  CompMat3 t = p2.t() * H * p1;
  return t[0];
}
comp_d crossprod(const Point& z1, const Point& z2,
                 const Point& w1, const Point& w2, const CompMat3& H)
{
  return (herm(w1, z1, H) * herm(w2, z2, H)) / 
         (herm(w2, z1, H) * herm(w1, z2, H));
}

double Ma(const CompMat3& A)
{
  return std::sqrt(std::abs(arma::trace(A) * arma::trace(A) - 4.0));
}

double comp_distance(const Point& p1, const Point& p2, const CompMat3& H)
{
  // this should have imaginary part 0
  // const comp_d t = (herm(p1, p2, H) * herm(p2, p1, H)) / 
  //                  (herm(p1, p1, H) * herm(p2, p2, H));
  const comp_d t = crossprod(p1, p2, p2, p1, H);

  return 2.0 * std::acosh(std::sqrt(std::real(t)));
}

Point find_null(const CompMat3& H)
{
  double lower_bound = -100;
  double upper_bound = 100; 
  std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
  std::random_device rd;
  std::mt19937 gen(rd());
  Point base(3);

  while (true)
  {
    base << comp_d(unif(gen), unif(gen)) << arma::endr
         << comp_d(unif(gen), unif(gen)) << arma::endr
         << comp_d(unif(gen), unif(gen)) << arma::endr;

    if (std::real(herm(base, base, H)) < 0)
      break;
  }
  return base;
}

bool IsWordMinimal(const unsigned int word,
                   const size_t n)
{
  // the identity is minimal!
  if (word == 0)
    return true;

  int w = word - 1;
  int prev = -1000; // big number
  size_t len = 2 * n;
  while(true)
  {
    const int e = w % len;
    if ((e % 2) == 0)
      if ((e + 1) == prev)
        return false;
    else
      if (e == (prev + 1))
        return false;

    prev = e;

    w /= len;
    if (w == 0)
      break;
    w -= 1;
  }
  return true;   
}

bool CheckWords2(const unsigned int upto,
                 const std::vector<CompMat3>& Mats,
                 const CompMat3& H,
                 bool& SeenParabolic)
{
  for (int word = 0; word < upto; ++word)
  {
    if (!IsWordMinimal(word, Mats.size()))
      continue;

    CompMat3 test = WordToMatrix(word, Mats);
    if (is_non_finite_elliptic(test, H))
    {
      // std::cout << std::endl;
      
      // arma::cx_vec eigval;
      // arma::cx_mat eigvec;

      // arma::eig_gen(eigval, eigvec, test);

      // std::cout << "evals = \n" << eigval <<std::endl;
      // std::cout << "cross = \n" << eigvec.t() * H * eigvec << std::endl;
      // std::cout << GetIsomClassStr(test, H) << std::endl;
      // std::cout << Order(test, H) << std::endl;
      return false;
    }

    const IsomClass iclass = GetIsomClass(test, H);

    if ((iclass == IsomClass::ParabolicPure) ||
        (iclass == IsomClass::ParabolicScrew))
    {
      SeenParabolic = true;
      continue;
    }
  }
  return true;  
}

void printShortWords(const unsigned int upto,
                     const std::vector<CompMat3>& Mats,
                     const std::vector<std::string>& Names,
                     const CompMat3& H,
                     const bool skip_loxo)
{
  unsigned int k = 1;
  for (unsigned int i = 0; i < upto; ++i)
  {
    if (!IsWordMinimal(i, Mats.size()))
      continue;

    CompMat3 temp = WordToMatrix(i, Mats);
    if ((GetIsomClass(temp, H) == IsomClass::Loxodromic) && skip_loxo)
      continue;

    std::cout << HumanWord(i, Names) << ": "
              << GetIsomClassStr(temp, H) << std::endl;
    ++k;
  }
}

Point getEllipticFixedPoint(const CompMat3& M, const CompMat3& H)
{
  arma::cx_vec eigval;
  arma::cx_mat eigvec;

  arma::eig_gen(eigval, eigvec, M);

  CompMat3 mHm = eigvec.t() * H * eigvec;

  const comp_d d1 = mHm(0,0);
  const comp_d d2 = mHm(1,1);
  const comp_d d3 = mHm(2,2);
  // make sure they are ALMOST real
  if ((std::imag(d1) > TOL) || 
      (std::imag(d2) > TOL) ||
      (std::imag(d3) > TOL))
    std::cout << "OH NO!" << std::endl;

  const double r1 = std::real(d1);
  const double r2 = std::real(d2);
  const double r3 = std::real(d3);

  // make sure exactly one is negative and two are positive
  const double pos = (r1 > TOL ? 1 : 0) +
                     (r2 > TOL ? 1 : 0) +
                     (r3 > TOL ? 1 : 0);

  const double neg = (r1 < TOL ? 1 : 0) +
                     (r2 < TOL ? 1 : 0) +
                     (r3 < TOL ? 1 : 0);

  if ((pos != 2) || (neg != 1))
    std::cout << "IT'S ALL GONE WRONG!" << std::endl;

  const size_t neg_index = (r1 < TOL ? 0 : (r2 < TOL ? 1 : 2));

  Point p(3);
  p << eigvec(0,neg_index) << arma::endr
    << eigvec(1,neg_index) << arma::endr
    << eigvec(2,neg_index) << arma::endr;

  // sanity checks
  // std::cout << herm(p, p, H) << std::endl;
  // std::cout << normalize(M * p) << std::endl;
  // std::cout << normalize(p) << std::endl;

  return normalize(p, 2);
}

Polar getLineReflectionPolar(const CompMat3& M, const CompMat3& H)
{
  //std::cout << "getLineReflectionPolar" << std::endl;
  arma::cx_vec eigval;
  arma::cx_mat eigvec;

  arma::eig_gen(eigval, eigvec, M);

  CompMat3 mHm = eigvec.t() * H * eigvec;

  const comp_d d0 = mHm(0,0);
  const comp_d d1 = mHm(1,1);
  const comp_d d2 = mHm(2,2);

  // make sure they are ALMOST real
  if ((std::imag(d0) > TOL) || 
      (std::imag(d1) > TOL) ||
      (std::imag(d2) > TOL))
    std::cout << "OH NO!" << std::endl;

  int odd_one_out = -1;
  arma::cx_dvec v_pos(3);

  if (std::abs(eigval[0] - eigval[1]) < TOL)
    odd_one_out = 2;
  else if (std::abs(eigval[1] - eigval[2]) < TOL)
    odd_one_out = 0;
  else if (std::abs(eigval[2] - eigval[0]) < TOL)
    odd_one_out = 1;
  
  v_pos << eigvec(0,odd_one_out) << arma::endr
        << eigvec(1,odd_one_out) << arma::endr
        << eigvec(2,odd_one_out) << arma::endr;

  return normalize(v_pos, 2);
}

void getLoxodromicFixed(const CompMat3& M, const CompMat3& H,
                        Point& p1, Point& p2)
{
  //std::cout << "getLineReflectionPolar" << std::endl;
  arma::cx_vec eigval;
  arma::cx_mat eigvec;

  arma::eig_gen(eigval, eigvec, M);

  CompMat3 mHm = eigvec.t() * H * eigvec;

  const comp_d d0 = mHm(0,0);
  const comp_d d1 = mHm(1,1);
  const comp_d d2 = mHm(2,2);

  int pos = -1;
  int null_1 = -1;
  int null_2 = -1;

  if (std::abs(std::abs(eigval[0]) - 1) < TOL)
  {
    pos = 0; null_1 = 1; null_2 = 2;
  }
  else if (std::abs(std::abs(eigval[1]) - 1) < TOL)
  {
    pos = 1; null_1 = 2; null_2 = 0;
  }
  else if (std::abs(std::abs(eigval[2]) - 1) < TOL)
  {
    pos = 2; null_1 = 0; null_2 = 1;
  }
  
  p1 << eigvec(0,null_1) << arma::endr
     << eigvec(1,null_1) << arma::endr
     << eigvec(2,null_1) << arma::endr;

  p2 << eigvec(0,null_2) << arma::endr
     << eigvec(1,null_2) << arma::endr
     << eigvec(2,null_2) << arma::endr;

  return;
}

void FindIsometryThatFixes(const Point& p,
                           const unsigned int upto,
                           const std::vector<CompMat3>& Mats,
                           const std::vector<std::string>& Names,
                           const CompMat3& H)
{
  for (unsigned int i = 0; i < upto; ++i)
  {
    if (!IsWordMinimal(i, Mats.size()))
      continue;
    CompMat3 temp = WordToMatrix(i, Mats);
    const Point Mp = temp * p;

    if (arma::approx_equal(normalize(Mp), normalize(p), "absdiff", TOL))
      std::cout << HumanWord(i, Names) << ": "
                << GetIsomClassStr(temp, H) << std::endl;
  }
}

CompMat3 WordToMatrix(const unsigned int word,
                      const std::vector<CompMat3>& Mats)
{
  int w = word - 1;
  size_t len = 2 * Mats.size();
  CompMat3 mat = arma::eye<CompMat3>(3,3);
  while(w >= 0)
  {
    const int e = w % len;
    if ((e % 2) == 0)
      mat *= Mats[e / 2];
    else
      mat *= arma::inv(Mats[e / 2]);
    w /= len;
    if (w == 0)
      break;
    w -= 1;
  }
  return mat;  
}

std::string HumanWord(const unsigned int word,
                      const std::vector<std::string>& Names)
{
  if (word == 0)
    return "Id";

  int w = word - 1;
  std::string s;
  size_t len = 2 * Names.size();
  while(w >= 0)
  {
    const int e = w % len;
    if ((e % 2) == 0)
      s.append(Names[e / 2]);
    else
    {
      s.append(Names[e / 2]);
      s.append("`");
    }
    w /= len;
    if (w == 0)
      break;
    w -= 1;
  }
  return s;
}

unsigned int StringToCode(const std::string word,
                          const std::vector<std::string>& Names)
{
  if (word == "Id")
    return 0;

  unsigned int w = 0;
  unsigned int pos = 0;

  size_t tot_names = Names.size();
  size_t len = word.length();

  int from = 0;
  while (from < len)
  {
    bool inverse = false;
    const std::string c = word.substr(from, 1);
    ++from;
    // check if next char is `
    const std::string inv = word.substr(from, 1);
    if (inv == "`")
    {
      inverse = true;
      ++from;
    }
    unsigned int i = 0;
    for (i = 0; i < Names.size(); ++i)
    {
      if (c == Names[i])
      {
        int bin = (2 * i) + (inverse ? 1 : 0) + 1;
        w += (std::pow(2 * tot_names, pos)) * bin;
        break;
      }
    }
    ++pos;
  }
  return w;
}


bool CheckWordsAlt(const std::vector<CompMat3>& Mats,
                   const CompMat3& H,
                   bool& SeenParabolic)
{
  for (size_t i = 0; i < Mats.size(); ++++i)
  {
    CompMat3 test = Mats[i];
    if (is_non_finite_elliptic(test, H))
      return false;

    const IsomClass iclass = GetIsomClass(test, H);

    if ((iclass == IsomClass::ParabolicPure) ||
        (iclass == IsomClass::ParabolicScrew))
    {
      SeenParabolic = true;
      continue;
    }
  }
  return true;
}

void getShortWords(const std::vector<CompMat3>& Generators,
                   std::vector<CompMat3>& Matrices,
                   const size_t upto)
{
  Matrices.reserve(upto);
  for (unsigned int i = 0; i <= upto; ++i)
  {
    if (!IsWordMinimal(i, Generators.size()))
       continue;

    const CompMat3 A = WordToMatrix(i, Generators);
    bool seen = false;
    int j = 0;
    for (j = 0; j < Matrices.size(); ++j)
    {
      const CompMat3 B = Matrices[j];
      seen = areEqual(A, B);
      if (seen)
        break;
    }

    if (!seen)
      Matrices.push_back(A);
  }
  return;
}
