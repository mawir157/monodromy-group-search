#include "word.h"

Word::Word(std::vector<Generator> gen_vec,
           CompMat3 matrix,
           CompMat3 H_matrix,
           IsomClass iso_class,
           comp_d trace,
           int order) : 
    m_gen_vec(gen_vec)
  , m_matrix(matrix)
  , m_H_matrix(H_matrix)
  , m_iso_class(iso_class)
  , m_trace(trace)
  , m_order(order) {}

Word::Word(CompMat3 H_matrix) :
    m_gen_vec(0)
  , m_matrix(arma::eye<CompMat3>(3,3))
  , m_H_matrix(H_matrix)
  , m_iso_class(IsomClass::Identity)
  , m_trace(comp_d(3.0, 0.0))
  , m_order(1) {}

std::string Word::as_string() const
{
  if (m_gen_vec.size() == 0) 
    return "Id";

  std::string word_str;
  for (size_t i = 0; i < m_gen_vec.size(); ++i)
  {
    Generator gen = m_gen_vec[i];
    switch (gen)
    {
      case Generator::R1: word_str.append("R1"); break;
      case Generator::R2: word_str.append("R2"); break;
      case Generator::R3: word_str.append("R3"); break;
      case Generator::E1: word_str.append("E1"); break;
      case Generator::E2: word_str.append("E2"); break;
      case Generator::E3: word_str.append("E3"); break;
      //
      case Generator::A:  word_str.append("A"); break;
      case Generator::Ai: word_str.append("A'"); break;
      case Generator::B:  word_str.append("B"); break;
      case Generator::Bi: word_str.append("B'"); break;
      case Generator::ID:
                 default: word_str.append("[]");
    }
  }
  return word_str;
}

bool Word::is_equal(const Word& wd, const bool debug) const
{
  if (m_iso_class != wd.get_isom_class())
    return false;

  if (m_order != wd.get_order())
    return false;

  if (!traceEqual(m_trace, wd.get_trace()))
    return false;

  if (!areEqual(m_matrix, wd.get_matrix()))
    return false;

  return true;
}

bool Word::is_equal_inverse(const Word& wd, const bool debug) const
{
  if (m_iso_class != wd.get_isom_class())
    return false;

  if (m_order != wd.get_order())
    return false;

  if (!areEqual(m_matrix, arma::inv(wd.get_matrix())))
    return false;

  return true;
}

std::string Word::str_isom_class() const
{
    switch (m_iso_class)
    {
      case IsomClass::Identity:   return "Identity";
      case IsomClass::Elliptic:   return "Elliptic";
      case IsomClass::ReflectionPoint: return "Reflection Point";
      case IsomClass::ReflectionLine: return "Reflection Line";
      case IsomClass::ParabolicPure:  return "Parabolic Pure";
      case IsomClass::ParabolicScrew:  return "Parabolic Screw";
      case IsomClass::Loxodromic: return "Loxodromic";
                         default: return "UNCLASSIFIED";
    }  
}

Generator Word::last_element() const
{
  if (m_iso_class == IsomClass::Identity)
    return Generator::ID;

  if (m_gen_vec.size() > 0)
    return m_gen_vec.back();
}

Generator Word::first_element() const
{
  if (m_iso_class == IsomClass::Identity)
    return Generator::ID;

  if (m_gen_vec.size() > 0)
    return m_gen_vec.front();
}

Word Word::invert() const
{
  if ((m_iso_class == IsomClass::Identity) || (m_gen_vec.size() == 0))
  {
    return Word(m_H_matrix);
  }

  std::vector<Generator> new_vector;
  new_vector.reserve(m_gen_vec.size());
  
  for (size_t i = 0; i < m_gen_vec.size(); ++i)
  {
    Generator gen_to_add = inverse(m_gen_vec[m_gen_vec.size() - 1 - i]);
    new_vector.push_back(gen_to_add);
  }

  CompMat3 new_matrix = arma::inv(m_matrix);
  IsomClass new_IsoClass = m_iso_class;
  comp_d new_trace = m_trace;
  int new_order = m_order;

  Word new_word(new_vector, new_matrix, m_H_matrix,
                new_IsoClass, new_trace, new_order);

  return new_word;
}

void Word::simplify_gen_vec()
{
  if (m_gen_vec.size() <= 1)
    return;

  bool changed = true;
  while (changed)
  {
    changed = false;    
    for (size_t i = 0; i < m_gen_vec.size()-1; ++i) 
    {
      if (m_gen_vec[i] == inverse(m_gen_vec[i+1]))
      {
        m_gen_vec.erase(m_gen_vec.begin() + i, m_gen_vec.begin() + i + 2);
        changed = true && (m_gen_vec.size() > 1);
        break;
      }    
    } 
  }
  return;
}

void Word::conjugate(const Word& P)
{

//    CompMat3               m_matrix;
//    std::vector<Generator> m_gen_vec;
//    IsomClass              m_iso_class;
//    comp_d                 m_trace;
//    int                    m_order;
  // create the new isometry matrix
  m_matrix = P.get_matrix() * m_matrix * arma::inv(P.get_matrix());
  // the isoclass, trace and order do not change
  std::vector<Generator> base_vector = m_gen_vec;
  std::vector<Generator> conj_vector = P.get_gen_vec();
  m_gen_vec.clear();

  m_gen_vec.reserve(base_vector.size() + 2 * P.get_gen_vec().size());

  for (size_t i = 0; i < conj_vector.size(); ++i)
    m_gen_vec.push_back(conj_vector[i]);

  for (size_t i = 0; i < base_vector.size(); ++i)
    m_gen_vec.push_back(base_vector[i]);

  for (size_t i = 0; i < conj_vector.size(); ++i)
  {
    Generator gen_to_add = inverse(conj_vector[conj_vector.size() - 1 - i]);
    m_gen_vec.push_back(gen_to_add);
  }

  simplify_gen_vec();
}

Word power(const Word& base_word, const unsigned int p)
{
  Word new_word = Word(base_word.get_H_matrix());
  for (unsigned int i=0; i < p; ++i)
    new_word = new_word * base_word;    

  return new_word;
}

Word conjugate(const Word& base_word, const Word& conj_word)
{
  CompMat3 new_mat = conj_word.get_matrix() * 
                     base_word.get_matrix() * 
                     arma::inv(conj_word.get_matrix());
  // the isoclass, trace and order do not change
  IsomClass i_class = base_word.get_isom_class();
  unsigned int order = base_word.get_order();
  comp_d trace = base_word.get_trace();

  std::vector<Generator> base_vector = base_word.get_gen_vec();
  std::vector<Generator> conj_vector = conj_word.get_gen_vec();
  std::vector<Generator> new_vec;

  new_vec.reserve(base_vector.size() + 2 * conj_word.get_gen_vec().size());

  for (size_t i = 0; i < conj_vector.size(); ++i)
    new_vec.push_back(conj_vector[i]);

  for (size_t i = 0; i < base_vector.size(); ++i)
    new_vec.push_back(base_vector[i]);

  for (size_t i = 0; i < conj_vector.size(); ++i)
  {
    Generator gen_to_add = inverse(conj_vector[conj_vector.size() - 1 - i]);
    new_vec.push_back(gen_to_add);
  }

  Word new_word(new_vec, new_mat, base_word.get_H_matrix(),
                i_class, trace, order);
  new_word.simplify_gen_vec();
  return new_word;
}

Word Word::operator*(const Word& wd) const
{
  CompMat3 new_mat = m_matrix * wd.get_matrix();

  IsomClass i_class = GetIsomClass(new_mat, m_H_matrix);

  unsigned int order = 0;
  if (i_class == IsomClass::Loxodromic ||
      i_class == IsomClass::ParabolicPure ||
      i_class == IsomClass::ParabolicScrew)
    order = -1;
  else if (i_class == IsomClass::Identity)
    order = 1;
  else
    order = Order(new_mat, m_H_matrix);

  comp_d trace = arma::trace(new_mat);

  std::vector<Generator> new_vec;

  new_vec.reserve(m_gen_vec.size() + wd.get_gen_vec().size());

  for (size_t i = 0; i < m_gen_vec.size(); ++i)
    new_vec.push_back(m_gen_vec[i]);

  for (size_t i = 0; i < wd.get_gen_vec().size(); ++i)
    new_vec.push_back(wd.get_gen_vec()[i]);

  Word new_word(new_vec, new_mat, m_H_matrix,
                i_class, trace, order);
  new_word.simplify_gen_vec();
  return new_word;
}

bool Word::operator<(const Word& wd) const
{
  if (m_gen_vec.size() < wd.word_length())
    return true;

  if (m_gen_vec.size() > wd.word_length())
    return false;
  
  std::vector<Generator> comp_vec = wd.get_gen_vec();
  // at this point the two words are the same length
  for (size_t i = 0; i < m_gen_vec.size(); ++i)
  {
    if (m_gen_vec[i] < comp_vec[i])
      return true;
    else if (m_gen_vec[i] > comp_vec[i])
      return false;
  }
  // if we get here then the words are equal
  return false;
}

bool Word::operator>(const Word& wd) const
{
  if (m_gen_vec.size() > wd.word_length())
    return true;

  if (m_gen_vec.size() < wd.word_length())
    return false;
  
  std::vector<Generator> comp_vec = wd.get_gen_vec();
  // at this point the two words are the same length
  for (size_t i = 0; i < m_gen_vec.size(); ++i)
  {
    if (m_gen_vec[i] > comp_vec[i])
      return true;
    else if (m_gen_vec[i] < comp_vec[i])
      return false;
  }
  // if we get here then the words are equal
  return false;
}

