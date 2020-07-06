#include "mode_file.h"

int mode_file(const std::string filepath)
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
    std::vector<Word> tg = gd.find_alt(new_Words, 2,
                                       groups[i].get_ref_order());
    
    std::cout << tg.size() << std::endl;

    if (gd.ok)
    {
      std::cout << gd.print() << std::endl;
      std::cout << "*  *  *  *  *  *  *  *  *  *  *  *  *  *" <<std::endl;
      gd.find_alt(tg, 2, tg[0].get_order());
      std::cout << gd.print() << std::endl;

      get_words_upto_n(8, tg, p_H, new_Words);
      summary(new_Words);
      FunDom fd(new_Words, H);
      const std::vector<Word> fd_words = fd.get_vertices(true);
    }
    else
    {
      std::cout << "FAILED TO FIND TRIANGLE GROUP DESCRIPTION" << std::endl;
    }

    // for (size_t i = 0; i < 10; ++i)
    // {
    //   Point base_point = find_neg(H);
    //   StochFunDom fd(base_point, H);
    //   fd.addPoints(new_Words);
    //   fd.stochastic_lattice(10000, 100.0, true);
    // }

    // std::cout << (jorgensen(new_Words) ? "*PASSED*" : "*FAILED*") << std::endl;
    std::cout << "******************************************" <<std::endl;
    std::cout << std::endl << std::endl;	
	}
	return 1;
}