#pragma once

#include "word.h"

std::vector<Word> get_words_upto_n(const unsigned int n,
	                               const std::vector<Word>& gens,
                                   const CompMat3& H,
	                               const bool verbose=true);
