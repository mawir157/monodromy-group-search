#include "mode_loop.h"

int mode_loop(const unsigned int denom_a, const unsigned int denom_b,
	            const std::string filelabel)
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
                  StochFunDom fd(base_point, H);
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

  return 1;
}