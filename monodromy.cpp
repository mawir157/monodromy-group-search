#include "includes.h"
#include "utils.h"
#include "parsedline.h"
#include "rational.h"
#include "matrixfns.h"
#include "groupdescription.h"
#include "fundom.h"
#include "geometry.h"
#include "polytope.h"
#include "jorgensen.h"
#include "getwordsupton.h"
#include "word.h"

#include <stdlib.h>

int main(int argc, char *argv[])
{
  RunMode mode = RunMode::unknown;
  std::string filepath = "NULL";
  unsigned int denom_a = 0;
  unsigned int denom_b = 0;

  if (argc < 1)
  {
    std::cout << "invalid options" << std::endl;
    return -1;
  }
  else
  {
    if (std::strcmp(argv[1], "-f") == 0)
    {
      std::cout << "f" << std::endl;
      mode = RunMode::file; // read groups from file;
      filepath = argv[2];
    }
    else if (std::strcmp(argv[1], "-l") == 0)
    {
      std::cout << "l" << std::endl;
      mode = RunMode::loop;
      denom_a = std::stoi(argv[2]);
      denom_b = std::stoi(argv[3]);
    }
    else
    {
      std::cout << "unrecognised options" << std::endl;
      return -1;
    }
  }

  if (mode == RunMode::file)
  {
    std::vector<ParsedLine> groups = parseFile(filepath);

    bool para = false;
    for (size_t i = 0; i < groups.size(); ++i)
    {
      auto [ A, B, H ] = groups[i].toMatrices();
      std::shared_ptr<const CompMat3> p_H = std::make_shared<const CompMat3>(H);

      mat_sig sig = get_mat_sig(H);
      std::vector<CompMat3> mats {A, B};

      std::vector<Generator> v_a {Generator::A};
      Word word_A(v_a, A, p_H, GetIsomClass(A,H), arma::trace(A), Order(A,H));
      Word word_Ai = word_A.invert();

      std::vector<Generator> v_b {Generator::B};
      Word word_B(v_b, B, p_H, GetIsomClass(B,H), arma::trace(B), Order(B,H));
      Word word_Bi = word_B.invert();

      std::vector<Word> gen_words {word_A, word_Ai, word_B, word_Bi};

      std::vector<Word> new_Words;
      bool good = get_words_upto_n(11, gen_words, p_H, new_Words);
      std::cout << (good ? "good" : "bad") << std::endl;

      std::cout << groups[i].print() << std::endl;
      std::cout << new_Words.size() << std::endl;

      std::cout << (CheckWords(new_Words, para) ? "" : " / NON-DISCRETE!");
      std::cout << (para ? " / PARABOLIC" : "") << std::endl;

      GroupDescription gd;
      gd.find_alt(new_Words, 2, groups[i].get_ref_order());
      if (gd.ok)
        std::cout << gd.print() << std::endl;


      // std::cout << (jorgensen(new_Words) ? "*PASSED*" : "*FAILED*") << std::endl;
      summary(new_Words);

      std::cout << std::endl << std::endl;
    }
  }

  if (mode == RunMode::loop)
  {

    const std::vector<Generator> v_a {Generator::A};
    const std::vector<Generator> v_b {Generator::B};

    for (size_t a_1 = 0; a_1 < denom_a; ++a_1)
    {
      for (size_t a_2 = a_1 + 1; a_2 < denom_a; ++a_2)
      {
        for (size_t a_3 = a_2 + 1; a_3 < denom_a; ++a_3)
        {
          const Triple ATriple(a_1,denom_a, a_2,denom_a, a_3,denom_a);
          const CompMat3 A = TripleToMatrix(ATriple);
          if (std::abs(arma::det(A)) < TOL)
            continue;

          for (size_t b_1 = 0; b_1 < denom_b; ++b_1)
          {
            for (size_t b_2 = b_1 + 1; b_2 < denom_b; ++b_2)
            {
              for (size_t b_3 = b_2 + 1; b_3 < denom_b; ++b_3)
              {
                const Triple BTriple(b_1,denom_b, b_2,denom_b, b_3,denom_b);
                const CompMat3 B = TripleToMatrix(BTriple);

                if (std::abs(arma::det(B)) < TOL)
                  continue;

                const CompMat3 H = SextupleToH(ATriple, BTriple);
                mat_sig sig = get_mat_sig(H);

                if (sig.match(2, 0, 1))
                  continue;

                std::shared_ptr<const CompMat3> p_H =
                                            std::make_shared<const CompMat3>(H);

                const Word word_A(v_a, A, p_H, GetIsomClass(A,H),
                                  arma::trace(A), Order(A,H));
                const Word word_Ai = word_A.invert();

                const Word word_B(v_b, B, p_H, GetIsomClass(B,H),
                                  arma::trace(B), Order(B,H));
                const Word word_Bi = word_B.invert();

                std::vector<Word> gen_Words {word_A, word_Ai, word_B, word_Bi};
                std::vector<Word> new_Words;
                bool ok = get_words_upto_n(9, gen_Words, p_H, new_Words, false);

                if (ok)
                {
                  std::cout << ATriple.asString() << " ";
                  std::cout << BTriple.asString() << std::endl;
                  std::cout << new_Words.size() << std::endl;
                  summary(new_Words);
                  GroupDescription gd;
                  gd.find_alt(new_Words, 2);
                  if (gd.ok)
                    std::cout << gd.print() << std::endl;

                  std::cout << std::endl;
                }
              }
            }
          }
        }
      }
    }
  }

  return 0;

////////////////////////////////////////////////////////////////////////////////
/*
  const Triple ATriple(5,14, 9,14, 7,14);
  const Triple BTriple(23,42, 24,42, 37,42);
  // const Triple ATriple(3,10, 7,10, 2,3);
  // const Triple BTriple(2,15, 5,15, 8,15);

  // const int alow = 15;
  // const int blow = 12;
  // srand(time(NULL));
  // const Triple ATriple(rand() % alow,alow,
  //                      rand() % alow,alow,
  //                      rand() % alow,alow);
  // const Triple BTriple(rand() % blow,blow,
  //                      rand() % blow,blow,
  //                      rand() % blow,blow);

  std::cout << "A = "<< ATriple.asString() << std::endl;
  std::cout << "B = "<< BTriple.asString() << std::endl;

  const CompMat3 A = TripleToMatrix(ATriple);
  const CompMat3 B = TripleToMatrix(BTriple);
  const CompMat3 H = SextupleToH(ATriple, BTriple);
  mat_sig sig = get_mat_sig(H);
  sig.print_sig();
  std::vector<CompMat3> mats {A, B};
  std::vector<std::string> names {"A", "B"};

  // find the fixed points of regular eliptic words
  std::vector<Point> fixed_points;
  std::vector<int> fixed_point_words;
  for (int word = 0; word < 100; ++word)
  {
    if (!IsWordMinimal(word, mats.size()))
      continue;

    CompMat3 temp = WordToMatrix(word, mats);
    if (GetIsomClass(temp) == IsomClass::Elliptic)
    {
      Polar pTemp = getEllipticFixedPoint(temp, H);
      fixed_points.push_back(pTemp);
      fixed_point_words.push_back(word);
    }
  }

  std::cout << fixed_points.size() <<std::endl;
  // find pairs of points contained in c-lines of reflections
  for (size_t i = 0; i < fixed_points.size(); ++i)
  {
    const Point p1 = fixed_points[i];
    const CompMat3 M1 = WordToMatrix(fixed_point_words[i], mats);
    for (size_t j = i + 1; j < fixed_points.size(); ++j)
    {
      const Point p2 = fixed_points[j];
      const CompMat3 M2 = WordToMatrix(fixed_point_words[j], mats);

      if (isPower(M1, M2) || isPower(M2, M1))
        continue;

      for (int word = 0; word < 1000; ++word)
      {
        if (!IsWordMinimal(word, mats.size()))
          continue;

        CompMat3 temp = WordToMatrix(word, mats);
        if (GetIsomClass(temp) == IsomClass::Reflection)
        {
          Polar polTemp = getLineReflectionPolar(temp, H);
          const double h1 = std::abs(herm(p1, polTemp, H));
          const double h2 = std::abs(herm(p2, polTemp, H));
          if ((h1 < TOL) && (h2 < TOL))
          {
            std::cout << HumanWord(fixed_point_words[i], names) << "~"
                      << HumanWord(fixed_point_words[j], names) << " C "
                      << HumanWord(word, names) << std::endl;
            break;
          }
        }
      }
    }
  }

  CompMat3 test = conj(B * arma::inv(A), B * B);
  std::cout << "A: " << GetIsomClassStr(A) << std::endl;
  std::cout << "B: " << GetIsomClassStr(B) << std::endl;
  CompMat3 BA = B * arma::inv(A);
  int ord = Order(BA);
  std::cout << "BA`: " << GetIsomClassStr(BA) << std::endl;
  std::cout << "test word: " << GetIsomClassStr(test) << std::endl;
  const Polar polBBBABB = getLineReflectionPolar(test, H);
  for (int word = 0; word < 1000; ++word)
  {
    if (!IsWordMinimal(word, mats.size()))
      continue;

    CompMat3 temp = WordToMatrix(word, mats);
    if (GetIsomClass(temp) == IsomClass::Elliptic)
    {
      Polar p1 = getEllipticFixedPoint(temp, H);
      const double h1 = std::abs(herm(p1, polBBBABB, H));
      if (h1 < TOL)
      {
        std::cout << HumanWord(word, names) << std::endl;
      }
    }
  }

  std::cout << A << std::endl;
  std::cout << test << std::endl;
  std::cout << "Braid(A, B) = " << Braid(A, B) << std::endl;
  std::cout << "Braid(A`, B) = " << Braid(arma::inv(A), B) << std::endl;
  std::cout << "Braid(A, test) = " << Braid(A, test) << std::endl;

  printShortWords(std::pow(4,9), mats, names);
  GroupDescription gd;
  gd.find(mats, names, 10, 1000000);
  std::cout << gd.print() << std::endl;

  std::cout << A << std::endl;
  std::cout << B * B * arma::inv(A) * B * arma::inv(B) * arma::inv(B) << std::endl;
  sig.print_sig();
  return 0;
*/
////////////////////////////////////////////////////////////////////////////////
/*
  std::string filename = GenerateFileName();
  std::cout << "Writing to " << filename << std::endl;

  std::ofstream output_file;
  output_file.open (filename);

  std::vector<Rational> rs = GetFirstRationals(100);
  std::cout << rs.size() << " ~ " << std::pow(1.0*rs.size(), 6.0) << std::endl;

  rs = GetFirstRationals(25);
  std::cout << rs.size() << " ~ " << std::pow(1.0*rs.size(), 6.0) << std::endl;

  rs = GetFirstRationals(26, false);
  std::cout << rs.size() << " ~ " << std::pow(1.0*rs.size(), 6.0) << std::endl;

  std::cout << std::endl;

  std::vector<Rational>::iterator it_a1;
  std::vector<Rational>::iterator it_a2;
  std::vector<Rational>::iterator it_a3;
  std::vector<Rational>::iterator it_b1;
  std::vector<Rational>::iterator it_b2;
  std::vector<Rational>::iterator it_b3;

  int count = 0;
  int seen = 0;
  const clock_t start = clock();

  for(it_a1 = rs.begin(); it_a1 != rs.end(); ++it_a1)
  {
    for(it_a2 = rs.begin(); it_a2 != rs.end(); ++it_a2)
    {
      if (*it_a2 <= *it_a1) { continue; }
      for(it_a3 = rs.begin(); it_a3 != rs.end(); ++it_a3)
      {
        if (*it_a3 <= *it_a2) { continue; }
        const Triple ATriple(*it_a1, *it_a2, *it_a3);
        const CompMat3 A = TripleToMatrix(ATriple);
        for(it_b1 = rs.begin(); it_b1 != rs.end(); ++it_b1)
        {
          for(it_b2 = rs.begin(); it_b2 != rs.end(); ++it_b2)
          {
            if (*it_b2 <= *it_b1) { continue; }
            for(it_b3 = rs.begin(); it_b3 != rs.end(); ++it_b3)
            {
              if (*it_b3 <= *it_b2) { continue; }

                ++count;
                if (count % 10000 == 0)
                {
                  const clock_t end = clock();
                  const double time = (double) (end-start) / CLOCKS_PER_SEC * 1000;

                  double tToEnd = (time / count)*std::pow(1.0*rs.size(), 6.0);
                  tToEnd /= 1000;
                  std::string  strToEnd = std::to_string(tToEnd);
                  strToEnd.append(" seconds");
                  if (tToEnd > 2 * 60)
                  {
                    tToEnd /= 60;
                    strToEnd = std::to_string(tToEnd);
                    strToEnd.append(" minutes");
                  }
                  if (tToEnd > 2 * 60)
                  {
                    tToEnd /= 60;
                    strToEnd = std::to_string(tToEnd);
                    strToEnd.append(" hours");
                  }
                  if (tToEnd > 2 * 24)
                  {
                    tToEnd /= 24;
                    strToEnd = std::to_string(tToEnd);
                    strToEnd.append(" days");
                  }
                  if (tToEnd > 2 * 365.25)
                  {
                    tToEnd /= 24;
                    strToEnd = std::to_string(tToEnd);
                    strToEnd.append(" years");
                  }
                  if (tToEnd > 2 * 100)
                  {
                    tToEnd /= 100;
                    strToEnd = std::to_string(tToEnd);
                    strToEnd.append(" centuries");
                  }

                  std::cout << "[" << count << " ("
                            << 100 * count / std::pow(1.0*rs.size(), 6.0)
                            << "%) ] " << "("
                            << (*it_a1).asString() << ", "
                            << (*it_a2).asString() << ", "
                            << (*it_a3).asString()
                            << " | "
                            << (*it_b1).asString() << ", "
                            << (*it_b2).asString() << ", "
                            << (*it_b3).asString()
                            << ") ["
                            << time / count << "ms per group, expected time to finish: "
                            << strToEnd
                            << "]" << std::endl;
                }

              const Triple BTriple(*it_b1, *it_b2, *it_b3);
              const CompMat3 B = TripleToMatrix(BTriple);

              std::vector<CompMat3> mats {A, B};
              std::vector<std::string> names {"A", "B"};
              const CompMat3 H = SextupleToH(ATriple, BTriple);
              mat_sig sig = get_mat_sig(H);
              if ((sig.m_pos != 2) || (sig.m_neg != 1))
                continue;
/////////////////////
              // std::cout << ATriple.asString()
              //           << " | "
              //           << BTriple.asString() << std::endl;
              // Point base = find_null(H);
              // FunDom fd(base, H);
              // unsigned int i = 0;
              // for (; i < std::pow(4, 9); ++i)
              // {
              //   CompMat3 temp = WordToMatrix(i, mats);
              //   fd.addPoint(temp);
              // }
              // fd.stochastic_lattice(1000, true);
/////////////////////
              // const CompMat3 gen = WordToMatrix(10, mats);
              // if (GetIsomClass(gen) != IsomClass::Reflection)
              //   continue;

              // int genOrder = Order(gen);

              // if ((genOrder < 2) || (genOrder > MAX_ORDER))
              //   continue;
              bool para = false;
//std::cout << "*";
              if (CheckWords2(std::pow(4, 6), mats, H, para))
              {
++seen;
//std::cout << std::endl;
                std::vector<CompMat3> words;
                getShortWords(mats, words, std::pow(4, 10));
//std::cout << "&" << std::endl;
                if (jorgensen(words, H))
                {
                  GroupDescription gd;
                  gd.find(mats, names, 10);
                  gd.even_better_find(mats, names, H, -1, 1000);
                  if (gd.ok)
                  {
                    std::cout << "("
                              << (*it_a1).asString() << ", "
                              << (*it_a2).asString() << ", "
                              << (*it_a3).asString()
                              << " | "
                              << (*it_b1).asString() << ", "
                              << (*it_b2).asString() << ", "
                              << (*it_b3).asString()
                              << ") | ";
                    std::cout << (1.0 * seen) / count << std::endl;
                    std::cout << gd.print() << std::endl;

//                    output_file << "("
//                                << (*it_a1).asString() << ", "
//                                << (*it_a2).asString() << ", "
//                                << (*it_a3).asString()
//                                << " | "
//                                << (*it_b1).asString() << ", "
//                                << (*it_b2).asString() << ", "
//                                << (*it_b3).asString()
//                                << ") | ";
//                    output_file << gd.print() << std::endl;
//                    std::vector<CompMat3> trimats {gd.R1, gd.R2, gd.R3};
//                    std::vector<std::string> trinames {"R1", "R2", "R3"};
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  output_file.close();

  return 0;
*/
}