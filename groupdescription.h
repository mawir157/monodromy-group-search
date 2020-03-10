#pragma once
#include "matrixfns.h"

class GroupDescription
{
  public:
    GroupDescription();
    CompMat3 R1;
    CompMat3 R2;
    CompMat3 R3;
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

    void find(const std::vector<CompMat3>& Mats,
              const std::vector<std::string>& Names,
              const CompMat3& H,
              const unsigned int seed,
              const unsigned int upto = MAX_WORD);
    void better_find(const std::vector<CompMat3>& Mats,
                     const std::vector<std::string>& Names,
                     const CompMat3& H,
                     const int ref_order,
                     const unsigned int upto = MAX_WORD);
    void even_better_find(const std::vector<CompMat3>& Mats,
                          const std::vector<std::string>& Names,
                          const CompMat3& H,
                          const int ref_order,
                          const unsigned int upto = MAX_WORD);
    void calc(const std::vector<CompMat3>& Mats,
              const std::vector<std::string>& Names,
              const CompMat3& H,
              const std::string w_1,
              const std::string w_2,
              const std::string w_3);
    void reset();
    std::string print();
};
