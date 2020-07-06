#include "includes.h"

// run modes
#include "mode_file.h"
#include "mode_loop.h"
#include "mode_matrix.h"

#include <stdlib.h>
#include <unistd.h> // used for cli options

/*
  -f read input from file
  -l max word length, default 10
  -a search a denom
  -b search b denom
  -v verbose mode - for debugging 
  -m read matrices from file
  [no arguement -> gtk]
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
    mode_file(filepath);
  }

  if (mode == RunMode::loop)
  {
    mode_loop(denom_a, denom_b, filelabel);
  }

  if (mode == RunMode::matrix)
  {
    mode_matrix(filepath);
  }

  if (mode == RunMode::verbose)
  {
    bool para = false;
    size_t index = 65;

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
