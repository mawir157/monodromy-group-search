#pragma once

#include "word.h"

bool get_words_upto_n(const unsigned int n,
	                    const std::vector<Word>& gens,
	                    std::shared_ptr<const CompMat3> H,
	                    std::vector<Word>& seen_words,
	                    const bool get_all = true);

bool CheckWords(const std::vector<Word>& Words,
                bool& SeenParabolic);

double summary(const std::vector<Word>& Words);