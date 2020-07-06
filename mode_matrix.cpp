#include "mode_matrix.h"

int mode_matrix(const std::string filepath)
{
	 auto [ H, mats] = parseMatrixFile(filepath);

	p_CompMat3 p_H = std::make_shared<const CompMat3>(H);

	std::vector<Generator> v_1 {Generator::R1};
	Word word_R1(v_1, mats[0], p_H, GetIsomClass(mats[0],H),
	            arma::trace(mats[0]), Order(mats[0],H));
	Word word_E1 = word_R1.invert();

	std::vector<Generator> v_2 {Generator::R2};
	Word word_R2(v_2, mats[1], p_H, GetIsomClass(mats[1],H),
	            arma::trace(mats[1]), Order(mats[1],H));
	Word word_E2 = word_R2.invert();

	std::vector<Generator> v_3 {Generator::R3};
	Word word_R3(v_3, mats[2], p_H, GetIsomClass(mats[2],H),
	            arma::trace(mats[2]), Order(mats[2],H));
	Word word_E3 = word_R3.invert();

	std::vector<Word> gen_words {word_R1, word_R2, word_R3,
	                             word_E1, word_E2, word_E3};

	const mat_sig sig = get_mat_sig(H);
	std::cout << sig.asString() << std::endl;

	std::vector<Word> new_Words;
	bool good = get_words_upto_n(15, gen_words, p_H, new_Words);
	std::cout << (good ? "DISCRETE!" : "NON-DISCRETE!") << std::endl;
	summary(new_Words);

	debugGroup(new_Words);

	GroupDescription gd(p_H);
	gd.find_alt(new_Words, 2);
	std::cout << gd.print() << std::endl;

	// for (size_t i = 0; i < 10; ++i)
	// {
	//   Point base_point = find_neg(H);
	//   StochFunDom sfd(base_point, H);
	//   sfd.addPoints(new_Words);
	//   sfd.stochastic_lattice(1000, 100.0, true);
	//   sfd.repel_lattice(100,100,true);
	// }

	FunDom fd(new_Words, H);
	const std::vector<Word> fd_words = fd.get_vertices(true);

	return 1;
}