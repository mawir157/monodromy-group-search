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
#include <unistd.h> // used for cli options

/*
  -f read input from file
  -l max word length, default 10
  -a search a denom
  -b search b denom
  -v verbose mode - for debugging 
*/


int main(int argc, char *argv[])
{
  RunMode mode = RunMode::unknown;
  std::string filepath = "NULL";
  unsigned int denom_a = 0;
  unsigned int denom_b = 0;

  unsigned int a1_1, a1_2, a2_1, a2_2, a3_1, a3_2;
  unsigned int b1_1, b1_2, b2_1, b2_2, b3_1, b3_2;

  std::string filelabel = "";

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
      filelabel.append("[");
      filelabel.append(argv[2]);
      filelabel.append("-");
      filelabel.append(argv[3]);
      filelabel.append("]");
    }
    else if (std::strcmp(argv[1], "-v") == 0)
    {
      std::cout << "v" << std::endl;
      mode = RunMode::verbose;
      filepath = argv[2];
    }
    else if (std::strcmp(argv[1], "-m") == 0)
    {
      std::cout << "m" << std::endl;
      mode = RunMode::matrix;
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

      p_CompMat3 p_H = std::make_shared<const CompMat3>(H);

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

      for (size_t i = 0; i < 10; ++i)
      {
        Point base_point = find_neg(H);
        FunDom fd(base_point, H);
        fd.addPoints(new_Words);
        fd.stochastic_lattice(10000, 100.0, true);
      }

      // std::cout << (jorgensen(new_Words) ? "*PASSED*" : "*FAILED*") << std::endl;

      std::cout << std::endl << std::endl;
    }
  }

  if (mode == RunMode::loop)
  {
    const std::vector<Generator> v_a {Generator::A};
    const std::vector<Generator> v_b {Generator::B};

    const std::string file_path = GenerateFileName("output", filelabel);
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

                p_CompMat3 p_H = std::make_shared<const CompMat3>(H);

                const Word word_A(v_a, A, p_H, GetIsomClass(A,H),
                                  arma::trace(A), Order(A,H));
                const Word word_Ai = word_A.invert();

                const Word word_B(v_b, B, p_H, GetIsomClass(B,H),
                                  arma::trace(B), Order(B,H));
                const Word word_Bi = word_B.invert();

                std::vector<Word> gen_Words {word_A, word_Ai, word_B, word_Bi};
                std::vector<Word> new_Words;
                bool ok = get_words_upto_n(10, gen_Words, p_H, new_Words, false);

                std::cout << (ok ? "DISCRETE!" : "NON-DISCRETE!") << std::endl;

                if (ok)
                {
                  std::cout << ATriple.asString() << " ";
                  std::cout << BTriple.asString() << std::endl;
                  summary(new_Words);
                  GroupDescription gd(p_H);
                  gd.find_alt(new_Words, 2, -1);
                  if (gd.ok)
                    std::cout << gd.print() << std::endl;
                  else
                    std::cout << "Cannot find a triangle group presentation"
                              << std::endl;

                  unsigned int lattice_count = 0;
                  for (size_t i = 0; i < 5; ++i)
                  {
                    Point base_point = find_neg(H);
                    FunDom fd(base_point, H);
                    fd.addPoints(new_Words);
                    if (fd.stochastic_lattice(10000, 100.0, true) < 100)
                      ++lattice_count;
                  }

                  output_file << ATriple.formal()
                              << BTriple.formal() << "*" << std::endl;
                  if (gd.ok)
                  {
                    output_file << "!";
                    output_file << gd.print() << std::endl;
                  }


                  if (lattice_count > 0)
                    output_file << "!" << "Lattice count = "
                                << lattice_count << std::endl;

                  std::cout << std::endl << std::endl;
                }
              }
            }
          }
        }
      }
    }
    output_file << "!" << "Loop finished without error" << std::endl;

    output_file.close();
  }

  if (mode == RunMode::matrix)
  {
    auto [ H, mats] = parseMatrixFile(filepath);
    std::cout << H << std::endl;

    for (size_t i = 0; i <  mats.size(); ++i)
    {
      std::cout << mats[i] << std::endl;
    }
  }

  if (mode == RunMode::verbose)
  {
    bool para = false;
    size_t index = 65;

//,8,30,9,30,23,30,
// (0, 1/10, 5/10) (3/10, 4/10, 8/10)
    const Triple ATriple(0,10, 1,10, 5,10);
    const CompMat3 A = TripleToMatrix(ATriple);

    // const Triple BTriple(8,30, 9,30, 23,30);
    const Triple BTriple(3,10, 4,10, 8,10);
    const CompMat3 B = TripleToMatrix(BTriple);

    const CompMat3 H = SextupleToH(ATriple, BTriple);
    const mat_sig sig = get_mat_sig(H);
    std::cout << sig.asString() << std::endl;

    // auto [ A, B, H ] = groups[index].toMatrices();
    p_CompMat3 p_H = std::make_shared<const CompMat3>(H);

    std::vector<CompMat3> mats {A, B};

    std::vector<Generator> v_a {Generator::A};
    Word word_A(v_a, A, p_H, GetIsomClass(A,H), arma::trace(A), Order(A,H));
    Word word_Ai = word_A.invert();

    std::vector<Generator> v_b {Generator::B};
    Word word_B(v_b, B, p_H, GetIsomClass(B,H), arma::trace(B), Order(B,H));
    Word word_Bi = word_B.invert();

    std::vector<Word> gen_words {word_A, word_Ai, word_B, word_Bi};

    std::vector<Word> new_Words;
    bool good = get_words_upto_n(10, gen_words, p_H, new_Words);
    std::cout << (good ? "DISCRETE!" : "NON-DISCRETE!") << std::endl;

    std::cout << ATriple.asString() << ", "
              << BTriple.asString() << std::endl;
    summary(new_Words);

    GroupDescription gd(p_H);
    gd.find_alt(new_Words, 2, -1);
    std::cout << gd.print() << std::endl;
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

    Word ell(p_H);
    // find a regular elliptic word
    for (size_t i = 0; i < new_Words.size(); ++i)
    {
      if (new_Words[i].get_isom_class() == IsomClass::Elliptic)
      {
        ell = new_Words[i];
        std::cout << new_Words[i].as_string() << std::endl;
        break;
      }
    }
    // get that word's fixed point
    Point p_ell = getEllipticFixedPoint(ell.get_matrix(), ell.get_H_matrix());

    std::vector<double> ds = get_spectrum(p_ell, new_Words);

    for (size_t i = 0; i < ds.size(); ++i)
    {
      std::cout << ds[i];
      if (i % 10 == 9)
        std::cout << std::endl;
      else
        std::cout << ", ";
    }

    // std::cout << (jorgensen(new_Words) ? "*PASSED*" : "*FAILED*") << std::endl;

    std::cout << std::endl << std::endl;
  }

  return 0;
}
