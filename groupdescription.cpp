#include "groupdescription.h"

GroupDescription::GroupDescription(std::shared_ptr<const CompMat3> H) :
    R1(H)
  , R2(H)
  , R3(H)
  , R1_ord(-163)
  , R2_ord(-163)
  , R3_ord(-163)
  , R123_ord(-163)
  , braid12(-163)
  , braid23(-163)
  , braid31(-163)
  , braid1323(-163)
  , braid2131(-163)
  , braid3212(-163)
  , ok(false)
{}

std::string GroupDescription::print()
{
  std::string s = "<";
  s.append(std::to_string(R1_ord));
  s.append("> (");
  s.append(std::to_string(braid23));
  s.append(",");
  s.append(std::to_string(braid31));
  s.append(",");
  s.append(std::to_string(braid12));
  s.append(";");
  s.append(std::to_string(braid1323));
  s.append(",");
  s.append(std::to_string(braid2131));
  s.append(",");
  s.append(std::to_string(braid3212));
  s.append("|");
  s.append(std::to_string(R123_ord));
  s.append(") <");
  s.append(R1.as_string());
  s.append(", ");
  s.append(R2.as_string());
  s.append(", ");
  s.append(R3.as_string());
  s.append(">");

  return s;
}

void GroupDescription::find_alt(const std::vector<Word>& words,
                                const int min_braid,
                                const int reflection_order)
{
  reset();
  std::vector<Word> reduced_words = remove_non_reflections(words);
  if (reduced_words.size() == 0)
  {
    std::cout << "No reflections found" << std::endl;
    return;
  }

  const CompMat3 H = reduced_words[0].get_H_matrix();

  std::vector<int> ref_ords;
  if (reflection_order == -1)
  {
    for (size_t i = 0; i < reduced_words.size(); ++i)
      ref_ords.push_back(reduced_words[i].get_order());

    std::sort(ref_ords.begin(), ref_ords.end());
    ref_ords.erase(std::unique(ref_ords.begin(), ref_ords.end()),
                   ref_ords.end());
  }
  else
  {
    ref_ords.push_back(reflection_order);
  }

  for (size_t i_ref = 0; i_ref < ref_ords.size(); ++i_ref)
  {
    const int ref_ord = ref_ords[i_ref];
    R1_ord = ref_ord; R2_ord = ref_ord; R3_ord = ref_ord;

    for (size_t i_1 = 0; i_1 < reduced_words.size(); ++i_1)
    {
      R1 = reduced_words[i_1];
      // make sure R1 has order ref_ord
      if (R1.get_order() != ref_ord)
        continue;

      for (size_t i_2 = i_1 + 1; i_2 < reduced_words.size(); ++i_2)
      {
        R2 = reduced_words[i_2];
        // make sure R2 has order ref_ord
        if (R2.get_order() != ref_ord)
          continue;

        const CompMat3 M1 = R1.get_matrix();
        const CompMat3 M2 = R2.get_matrix();

        // make sure they are not powers of each other
        if (isPower(M1, M2) || isPower(M2, M1))
            continue;

        // make sure they braid
        braid12 = Braid(M1, M2);
        if (braid12 < 3)
          continue;

        for (size_t i_3 = i_2 + 1; i_3 < reduced_words.size(); ++i_3)
        {
          R3 = reduced_words[i_3];
          // make sure R3 has order ref_ord
          if (R3.get_order() != ref_ord)
            continue;

          const CompMat3 M3 = R3.get_matrix();
          if (isPower(M1, M3) || isPower(M3, M1) ||
              isPower(M3, M2) || isPower(M2, M3))
            continue;

          braid31 = Braid(M1, M3);
          if (braid31 < 3)
            continue;

          braid23 = Braid(M3, M2);
          if (braid23 < 3)
            continue;

          const CompMat3 M323 = conj(arma::inv(M3), M2); // R3^-1*R2*R3
          braid1323 = Braid(M1, M323);
          if (braid1323 < min_braid)
            continue;

          const CompMat3 M131 = conj(arma::inv(M1), M3); // R1^-1*R3*R1
          braid2131 = Braid(M2, M131);
          if (braid2131 < min_braid)
            continue;

          const CompMat3 M212 = conj(arma::inv(M2), M1); // R2^-1*R1*R2
          braid3212 = Braid(M3, M212);
          if (braid3212 < min_braid)
            continue;

          const CompMat3 M123 = M1 * M2 * M3;
          R123_ord = Order(M123, H);
          if (R123_ord < 0)
            continue;

          // we've got here, so we have found a group.
          ok = true;
          // word1 = R1.as_string();
          // word2 = R2.as_string();
          // word3 = R3.as_string();
          return;
        }
      }
    }
  }
  // if we've got here, we have failed to find a group presentation
  // so reset and bail out
  reset();
  return;
}

std::vector<Word> GroupDescription::remove_non_reflections(const std::vector<Word>& words,
                                                           const bool allow_ref_points)
{
  std::vector<Word> v;
  for (size_t i = 0; i < words.size(); ++i)
  {
    if ((words[i].get_isom_class() == IsomClass::ReflectionLine) ||
        (allow_ref_points &&
         (words[i].get_isom_class() == IsomClass::ReflectionPoint)))
    {
      v.push_back(words[i]);
    }
  }
  return v;
}

void GroupDescription::reset()
{
  R1 = Word(R1.p_H_matrix());
  R2 = Word(R1.p_H_matrix());
  R2 = Word(R1.p_H_matrix());
  R1_ord = -163;
  R2_ord = -163;
  R3_ord = -163;
  R123_ord = -163;
  braid12 = -163;
  braid23 = -163;
  braid31 = -163;
  braid1323 = -163;
  braid2131 = -163;
  braid3212 = -163;
  ok = false;
}
