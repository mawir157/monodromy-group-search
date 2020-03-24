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

  unsigned int a1_1, a1_2, a2_1, a2_2, a3_1, a3_2;
  unsigned int b1_1, b1_2, b2_1, b2_2, b3_1, b3_2;

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
      std::cout << "searching all groups with of the form "
                << "*/" << denom_a << ", "
                << "*/" << denom_b << std::endl;
    }
    else if (std::strcmp(argv[1], "-v") == 0)
    {
      std::cout << "v" << std::endl;
      mode = RunMode::verbose;
      filepath = argv[2];
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

    for (size_t i = 0; i < groups.size(); ++i)
    {
      auto [ A, B, H ] = groups[i].toMatrices();
      mat_sig sig = get_mat_sig(H);

      if (!sig.match(2, 0, 1))
        continue;

      std::shared_ptr<const CompMat3> p_H = std::make_shared<const CompMat3>(H);

      std::vector<CompMat3> mats {A, B};

      std::vector<Generator> v_a {Generator::A};
      Word word_A(v_a, A, p_H, GetIsomClass(A,H), arma::trace(A), Order(A,H));
      Word word_Ai = word_A.invert();

      std::vector<Generator> v_b {Generator::B};
      Word word_B(v_b, B, p_H, GetIsomClass(B,H), arma::trace(B), Order(B,H));
      Word word_Bi = word_B.invert();

      std::vector<Word> gen_words {word_A, word_Ai, word_B, word_Bi};

      std::vector<Word> new_Words;
      std::cout << groups[i].print() << std::endl;
      bool good = get_words_upto_n(11, gen_words, p_H, new_Words);
      std::cout << (good ? "DISCRETE!" : "NON-DISCRETE!") << std::endl;
      summary(new_Words);

      GroupDescription gd(p_H);
      gd.find_alt(new_Words, 2, groups[i].get_ref_order());
      if (gd.ok)
        std::cout << gd.print() << std::endl;

      // std::cout << (jorgensen(new_Words) ? "*PASSED*" : "*FAILED*") << std::endl;

      std::cout << std::endl << std::endl;
    }
  }

  if (mode == RunMode::loop)
  {

    const std::vector<Generator> v_a {Generator::A};
    const std::vector<Generator> v_b {Generator::B};

    const std::string file_path = GenerateFileName();
    std::ofstream output_file;
    output_file.open(file_path);
    std::cout << "Saving to " << file_path << std::endl;

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

                CompMat3 H;
                try
                {
                  H = SextupleToH(ATriple, BTriple);
                }
                catch (const std::exception& e)
                {
                  std::cout << "could not build H" << std::endl;
                  continue;
                }

                // const CompMat3 H = SextupleToH(ATriple, BTriple);
                mat_sig sig = get_mat_sig(H);

                if (!sig.match(2, 0, 1))
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
                  GroupDescription gd(p_H);
                  gd.find_alt(new_Words, 2);
                  if (gd.ok)
                    std::cout << gd.print() << std::endl;

                  std::cout << std::endl;

                  output_file << ATriple.formal()
                              << BTriple.formal() << "*" << std::endl;
                  if (gd.ok)
                  {
                    output_file << "!";
                    output_file << gd.print() << std::endl;
                  }
                }
              }
            }
          }
        }
      }
    }

    output_file.close();
  }

  if (mode == RunMode::verbose)
  {
    bool para = false;
    std::vector<ParsedLine> groups = parseFile(filepath);

    auto [ A, B, H ] = groups[0].toMatrices();
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
    bool good = get_words_upto_n(12, gen_words, p_H, new_Words);
    std::cout << (good ? "DISCRETE!" : "NON-DISCRETE!") << std::endl;

    std::cout << groups[0].print() << std::endl;
    std::cout << new_Words.size() << std::endl;

    GroupDescription gd(p_H);
    gd.find_alt(new_Words, 2, -1);
    if (gd.ok)
    {
      std::cout << gd.print() << std::endl;
      std::cout << (gd.R1 * gd.R2).get_matrix() << std::endl;
      std::cout << (gd.R2 * gd.R1).get_matrix() << std::endl;
      std::cout << "==================================" << std::endl;
      std::cout << Braid(gd.R1.get_matrix(), gd.R2.get_matrix()) << std::endl;
      std::cout << "==================================" << std::endl;
      std::cout << (gd.R1 * gd.R1 * gd.R1 * gd.R1).get_matrix() << std::endl;
    }


    // std::cout << (jorgensen(new_Words) ? "*PASSED*" : "*FAILED*") << std::endl;
    summary(new_Words);

    std::cout << std::endl << std::endl;
  }

  return 0;
}
