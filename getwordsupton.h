#pragma once

#include "word.h"

std::vector<Word> get_words_upto_n(const unsigned int n,
	                               const std::vector<Word>& gens,
                                   const CompMat3& H,
	                               const bool verbose=true);

bool CheckWords(const std::vector<Word>& Words,
                bool& SeenParabolic);

double summary(const std::vector<Word>& Words);