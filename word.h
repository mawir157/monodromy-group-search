#pragma once

#include "matrixfns.h"

class Word
{
  public:
    // empty constructor creates identity word
    Word(std::shared_ptr<const CompMat3> H_matrix);
    // full constructor (probably never used, except for generators)
    Word(std::vector<Generator> gen_vec,
         CompMat3 matrix,
         std::shared_ptr<const CompMat3> H_matrix,
         IsomClass iso_class,
         comp_d trace,
         int order);
    // construct from vector
    explicit Word(std::vector<Generator> gen_vec);

    //
    std::string as_string() const;

    // check if it is equal to another Word
    bool is_equal(const Word& wd, const bool debug=false) const;
    // check if a word is equal to the inverse of another word
    bool is_equal_inverse(const Word& wd, const bool debug=false) const;

    // get the last generator of the word
    Generator last_element() const;
    // get the first generator of the word
    Generator first_element() const;
    // get the length of the word
    size_t word_length() const { return m_gen_vec.size(); }
    Word invert() const;
    void simplify_gen_vec();
    void overwrite_gen_vec(const std::vector<Generator>& v) { m_gen_vec = v;}

    // inplace conjugation
    void conjugate(const Word& P);

    // getters
    int get_order() const { return m_order; }
    IsomClass get_isom_class() const { return m_iso_class; }
    std::string str_isom_class() const;
    CompMat3 get_matrix() const { return m_matrix; }
    const CompMat3 get_H_matrix() const { return *(m_H_matrix.get()); }
    p_CompMat3 p_H_matrix() const { return m_H_matrix; }
    comp_d get_trace() const { return m_trace; }
    std::vector<Generator> get_gen_vec() const { return m_gen_vec; }
    bool is_non_finite_elliptic() const;

    // overload operators here
    Word operator*(const Word& wd) const;
    bool operator<(const Word& wd) const;
    bool operator>(const Word& wd) const;

    // image of point
    Point image(const Point p) const;

  private:
    CompMat3               m_matrix;
    p_CompMat3             m_H_matrix;
    std::vector<Generator> m_gen_vec;
    IsomClass              m_iso_class;
    comp_d                 m_trace;
    int                    m_order;
};

Word conjugate(const Word& base_word, const Word& conj_word);
Word power(const Word& base_word, const unsigned int p);

void FindIsometryThatFixes(const Point& p,
                           const std::vector<Word>& words);
