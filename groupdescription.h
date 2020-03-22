#pragma once

#include "word.h"

class GroupDescription
{
  public:
    GroupDescription();
    CompMat3 M1;
    CompMat3 M2;
    CompMat3 M3;
    int R1_ord;
    int R2_ord;
    int R3_ord;
    int R123_ord;
    std::string word1;
    std::string word2;
    std::string word3;
    int braid12;
    int braid23;
    int braid31;
    int braid1323;
    int braid2131;
    int braid3212;
    bool ok;

    void find_alt(const std::vector<Word>& words,
                  const int min_braid,
                  const int reflection_order = -1);
    void reset();
    std::string print();
  private:
    std::vector<Word> remove_non_reflections(const std::vector<Word>& words,
                                             const bool allow_ref_points = false);
};
