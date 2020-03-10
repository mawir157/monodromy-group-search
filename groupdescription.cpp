#include "groupdescription.h"

GroupDescription::GroupDescription() :
    R1(arma::eye<CompMat3>(3,3))
  , R2(arma::eye<CompMat3>(3,3))
  , R3(arma::eye<CompMat3>(3,3))
  , R1_ord(-163)
  , R2_ord(-163)
  , R3_ord(-163)
  , R123_ord(-163)
  , word1("UNSET")
  , word2("UNSET")
  , word3("UNSET")
  , braid12(-163)
  , braid23(-163)
  , braid31(-163)
  , braid1323(-163)
  , braid2131(-163)
  , braid3212(-163)
  , ok(false)
{}

void GroupDescription::find(const std::vector<CompMat3>& Mats,
                            const std::vector<std::string>& Names,
                            const CompMat3& H,
                            const unsigned int seed,
                            const unsigned int upto)
{
  R1 = WordToMatrix(seed, Mats);//B * arma::inv(A);
  R1_ord = Order(R1, H);
  word1 = HumanWord(seed, Names);//"BA`";

  for (int word_i2 = seed + 1; word_i2 < seed + upto; ++word_i2)
  {
    if (!IsWordMinimal(word_i2, Mats.size()))
      continue;

    R2 = WordToMatrix(word_i2, Mats);

    if (GetIsomClass(R2, H) != IsomClass::ReflectionLine)
      continue;

    if (isPower(R1, R2) || isPower(R2, R1))
      continue;

    R2_ord = Order(R2, H);
    word2 = HumanWord(word_i2, Names);
    if (R2_ord != R1_ord)
      continue;

    braid12 = Braid(R1, R2);
    if ((braid12) < 1)
      continue;

    // at this point we have R1, R2 reflections of the same order that braid 
    for (int word_i3 = word_i2 + 1; word_i3 < word_i2 + upto; ++word_i3)
    {
      if (!IsWordMinimal(word_i3, Mats.size()))
        continue;

      R3 = WordToMatrix(word_i3, Mats);

      if (GetIsomClass(R3, H) != IsomClass::ReflectionLine)
        continue;

      if (isPower(R1, R3) || isPower(R3, R1) ||
          isPower(R3, R2) || isPower(R2, R3))
        continue;

      R3_ord = Order(R3, H);
      word3 = HumanWord(word_i3, Names);
      if (R3_ord != R1_ord)
        continue;

      braid31 = Braid(R1, R3);  
      if (braid31 < 1)
        continue;

      braid23 = Braid(R3, R2);
      if (braid23 < 1)
        continue;

      const CompMat3 R323 = conj(arma::inv(R3), R2); // R3^-1*R2*R3
      braid1323 = Braid(R1, R323);
      if (braid1323 <= 1)
        continue;

      const CompMat3 R131 = conj(arma::inv(R1), R3); // R1^-1*R3*R1
      braid2131 = Braid(R2, R131);
      if (braid2131 <= 1)
        continue;

      const CompMat3 R212 = conj(arma::inv(R2), R1); // R2^-1*R1*R2
      braid3212 = Braid(R3, R212);
      if (braid3212 <= 1)
        continue;

      CompMat3 R123 = R1 * R2 * R3;
      R123_ord = Order(R123, H);
      if (R123_ord < 0)
        continue;

      ok = true;

      return;
    }
  }

  // failed to find reset
  reset();
}

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
  s.append(word1);
  s.append(", ");
  s.append(word2);
  s.append(", ");
  s.append(word3);
  s.append(">");
/*
  std::vector<std::string> names {"A", "B"};
  s.append(" <");
  s.append(std::to_string(StringToCode(word1, names)));
  s.append(", ");
  s.append(std::to_string(StringToCode(word2, names)));
  s.append(", ");
  s.append(std::to_string(StringToCode(word3, names)));
  s.append(">");
*/
  return s;
}

// we do not provide a seed word
void GroupDescription::better_find(const std::vector<CompMat3>& Mats,
                                   const std::vector<std::string>& Names,
                                   const CompMat3& H,
                                   const int ref_order,
                                   const unsigned int upto) 
{
  for (int word_i1 = 1; word_i1 < upto; ++word_i1)
  {
    if (!IsWordMinimal(word_i1, Mats.size()))
      continue;

    R1 = WordToMatrix(word_i1, Mats);

    const IsomClass isom_1 = GetIsomClass(R1, H);

    if ((isom_1 != IsomClass::ReflectionLine) &&
        (isom_1 != IsomClass::ReflectionPoint))
      continue;

    R1_ord = Order(R1, H);

    if (ref_order > 0)
      if (R1_ord != ref_order)
        continue;
    else
      if (R1_ord <= 0)
        continue;

    word1 = HumanWord(word_i1, Names);
//std::cout << "\t" << HumanWord(word_i1, Names) << std::endl;
    for (int word_i2 = word_i1 + 1; word_i2 < upto; ++word_i2)
    {
      if (!IsWordMinimal(word_i2, Mats.size()))
        continue;

      R2 = WordToMatrix(word_i2, Mats);
      const IsomClass isom_2 = GetIsomClass(R2, H);
      if (isom_2 != isom_1)
        continue;

      R2_ord = Order(R2, H);
      if (R2_ord != R1_ord)
        continue;

      if (isPower(R1, R2) || isPower(R2, R1))
        continue;

      braid12 = Braid(R1, R2);
      if (braid12 < 1)
        continue;

      word2 = HumanWord(word_i2, Names);
//std::cout << "\t\t" << HumanWord(word_i1, Names) << "-" << HumanWord(word_i2, Names) << std::endl;
      // at this point we have R1, R2 reflections of the same order that braid 
      for (int word_i3 = word_i2 + 1; word_i3 < upto; ++word_i3)
      {
// std::cout << "\t\t\t" << word_i3 << std::endl;
        if (!IsWordMinimal(word_i3, Mats.size()))
          continue;

        R3 = WordToMatrix(word_i3, Mats);
        const IsomClass isom_3 = GetIsomClass(R3, H);
        if (isom_3 != isom_1)
          continue;

        R3_ord = Order(R3, H);
        if (R3_ord != R1_ord)
          continue;

        if (isPower(R1, R3) || isPower(R3, R1) ||
            isPower(R3, R2) || isPower(R2, R3))
          continue;

        braid31 = Braid(R1, R3);  
        if (braid31 < 1)
          continue;

        braid23 = Braid(R3, R2);
        if (braid23 < 1)
          continue;

        const CompMat3 R323 = conj(arma::inv(R3), R2); // R3^-1*R2*R3
        braid1323 = Braid(R1, R323);
        if (braid1323 < 2)
          continue;

        const CompMat3 R131 = conj(arma::inv(R1), R3); // R1^-1*R3*R1
        braid2131 = Braid(R2, R131);
        if (braid2131 < 2)
          continue;

        const CompMat3 R212 = conj(arma::inv(R2), R1); // R2^-1*R1*R2
        braid3212 = Braid(R3, R212);
        if (braid3212 < 2)
          continue;

        CompMat3 R123 = R1 * R2 * R3;
        R123_ord = Order(R123, H);
        if (R123_ord < 0)
          continue;

        word3 = HumanWord(word_i3, Names);

        //std::cout << "\t" << braid1323 << "\t" << braid2131 << "\t" << braid3212 << std::endl;
        //std::cout << "\t" << (braid1323 < 2) << "\t" << (braid2131 < 2) << "\t" << (braid3212 < 2) << std::endl;

        //std::cout << "\t" << Braid(R1, R323) << "\t" << Braid(R2, R131) << "\t" << Braid(R3, R212) << std::endl;

        ok = true;

        return;
      }
    }
  }

  // failed to find reset
  reset();
  return;
}

void GroupDescription::even_better_find(const std::vector<CompMat3>& Mats,
                                        const std::vector<std::string>& Names,
                                        const CompMat3& H,
                                        const int ref_order,
                                        const unsigned int upto)
{
  // hash the matrix, order and isometry class as we see them so we don't have
  // to repeatedly recalculate them.
  std::map<int, CompMat3> matrices;
  std::map<int, IsomClass> iso_classes;
  std::map<int, int> orders;


  for (int cantor = 6; cantor < upto; ++cantor)
  {
    for (int word_i1 = 1; word_i1 < cantor; ++word_i1)
    {

      if (!IsWordMinimal(word_i1, Mats.size()))
        continue;

      if (!matrices.count(word_i1) > 0)
        matrices[word_i1] = WordToMatrix(word_i1, Mats);

      //R1 = WordToMatrix(word_i1, Mats);
      R1 = matrices[word_i1];

      if (!iso_classes.count(word_i1) > 0)
        iso_classes[word_i1] = GetIsomClass(R1, H);

      //const IsomClass isom_1 = GetIsomClass(R1, H);
      const IsomClass isom_1 = iso_classes[word_i1];

      if ((isom_1 != IsomClass::ReflectionLine) &&
          (isom_1 != IsomClass::ReflectionPoint))
        continue;

      if (!orders.count(word_i1) > 0)
        orders[word_i1] = Order(R1, H);

      //R1_ord = Order(R1, H);
      R1_ord = orders[word_i1];

      if (ref_order > 0)
        if (R1_ord != ref_order)
          continue;
      else
        if (R1_ord <= 0)
          continue;

      word1 = HumanWord(word_i1, Names);

      for (int word_i2 = 1; word_i1 + word_i2 < cantor; ++word_i2)
      {

        if (word_i1 == word_i2)
          continue;

          if (!IsWordMinimal(word_i2, Mats.size()))
            continue;

          if (!matrices.count(word_i2) > 0)
            matrices[word_i2] = WordToMatrix(word_i2, Mats);

          //R2 = WordToMatrix(word_i1, Mats);
          R2 = matrices[word_i2];

          if (!iso_classes.count(word_i2) > 0)
            iso_classes[word_i2] = GetIsomClass(R2, H);

          //const IsomClass isom_2 = GetIsomClass(R2, H);
          const IsomClass isom_2 = iso_classes[word_i2];
          if (isom_2 != isom_1)
            continue;

          if (!orders.count(word_i2) > 0)
            orders[word_i2] = Order(R2, H);

          //R2_ord = Order(R2, H);
          R2_ord = orders[word_i2];
          if (R2_ord != R1_ord)
            continue;

          if (isPower(R1, R2) || isPower(R2, R1))
            continue;

          braid12 = Braid(R1, R2);
          if (braid12 < 1)
            continue;

          word2 = HumanWord(word_i2, Names);

          int word_i3 = cantor - word_i2 - word_i1;

          if ((word_i3 == word_i2) || (word_i3 == word_i1))
            continue;

          if (!IsWordMinimal(word_i3, Mats.size()))
            continue;

          if (!matrices.count(word_i3) > 0)
            matrices[word_i3] = WordToMatrix(word_i3, Mats);

          //R3 = WordToMatrix(word_i3, Mats);
          R3 = matrices[word_i3];
          if (!iso_classes.count(word_i3) > 0)
            iso_classes[word_i3] = GetIsomClass(R3, H);

          //const IsomClass isom_3 = GetIsomClass(R3, H);
          const IsomClass isom_3 = iso_classes[word_i3];
          if (isom_3 != isom_1)
            continue;

          if (!orders.count(word_i3) > 0)
            orders[word_i3] = Order(R3, H);

          //R3_ord = Order(R3, H);
          R3_ord = orders[word_i3];
          if (R3_ord != R1_ord)
            continue;

          if (isPower(R1, R3) || isPower(R3, R1) ||
              isPower(R3, R2) || isPower(R2, R3))
            continue;

          braid31 = Braid(R1, R3);  
          if (braid31 < 3)
            continue;

          braid23 = Braid(R3, R2);
          if (braid23 < 3)
            continue;

          const CompMat3 R323 = conj(arma::inv(R3), R2); // R3^-1*R2*R3
          braid1323 = Braid(R1, R323);
          if (braid1323 < 3)
            continue;

          const CompMat3 R131 = conj(arma::inv(R1), R3); // R1^-1*R3*R1
          braid2131 = Braid(R2, R131);
          if (braid2131 < 3)
            continue;

          const CompMat3 R212 = conj(arma::inv(R2), R1); // R2^-1*R1*R2
          braid3212 = Braid(R3, R212);
          if (braid3212 < 2)
            continue;

          CompMat3 R123 = R1 * R2 * R3;
          R123_ord = Order(R123, H);
          if (R123_ord < 0)
            continue;

          word3 = HumanWord(word_i3, Names);

          return;   
      }
    }
  }

  // failed to find reset
  reset();
  return;
}


void GroupDescription::calc(const std::vector<CompMat3>& Mats,
                            const std::vector<std::string>& Names,
                            const CompMat3& H,
                            const std::string w_1,
                            const std::string w_2,
                            const std::string w_3)
{
  const int code1 = StringToCode(w_1, Names);
  const int code2 = StringToCode(w_2, Names);
  const int code3 = StringToCode(w_3, Names);

  R1 = WordToMatrix(code1, Mats);
  R2 = WordToMatrix(code2, Mats);
  R3 = WordToMatrix(code3, Mats);
  R1_ord = Order(R1, H);
  R2_ord = Order(R2, H);
  R3_ord = Order(R3, H);
  R123_ord = Order(R1 * R2 * R3, H);
  word1 = HumanWord(code1, Names);
  word2 = HumanWord(code2, Names);
  word3 = HumanWord(code3, Names);
  const CompMat3 R323 = conj(arma::inv(R3), R2); // R3^-1*R2*R3
  braid1323 = Braid(R1, R323);
  const CompMat3 R131 = conj(arma::inv(R1), R3); // R1^-1*R3*R1
  braid2131 = Braid(R2, R131);
  const CompMat3 R212 = conj(arma::inv(R2), R1); // R2^-1*R1*R2
  braid3212 = Braid(R3, R212);
  braid12 = Braid(R1, R2);
  braid31 = Braid(R1, R3); 
  braid23 = Braid(R3, R2);
  ok = true;
}

void GroupDescription::reset()
{
  R1 = arma::eye<CompMat3>(3,3);
  R2 = arma::eye<CompMat3>(3,3);
  R3 = arma::eye<CompMat3>(3,3);
  R1_ord = -163;
  R2_ord = -163;
  R3_ord = -163;
  R123_ord = -163;
  word1 = "*";
  word2 = "*";
  word3 = "*";
  braid12 = -163;
  braid23 = -163;
  braid31 = -163;
  braid1323 = -163;
  braid2131 = -163;
  braid3212 = -163;
  ok = false;
}