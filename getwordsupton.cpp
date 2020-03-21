#include "getwordsupton.h"

std::vector<Word> get_words_upto_n(const unsigned int n,
                                   const std::vector<Word>& gens,
                                   const CompMat3& H,
                                   const bool verbose)
{
  size_t i = 0;
  std::vector<Word> seen_words;
  seen_words.reserve(1000000);
  seen_words.push_back(Word(H));

  while (i < n)
  {
    std::vector<Word> seen_words_n_minus_1;
    seen_words_n_minus_1.reserve(1000000);
    for (size_t j = 0; j < seen_words.size(); ++j)
    {
      const Word test_word = seen_words[j];
      if (test_word.word_length() == i)
        seen_words_n_minus_1.push_back(test_word);
    }

    for (size_t j = 0; j < seen_words_n_minus_1.size(); ++j)
    {
      const Word base_word = seen_words_n_minus_1[j];
      const Generator base_last = base_word.last_element();

      for (size_t k = 0; k < gens.size(); ++k)
      {
        const Word new_gen = gens[k];
        // if the final element is the inverse of the generator we're about to
        // do nothing as they will cancel
        if (inverse(base_last) == new_gen.last_element())
          continue;

        // get a new word by multiplying the generator to the base word
        const Word new_word = base_word * new_gen;

        // check if we've seen this before!
        bool seen = false;
        for (unsigned int w = 0; w < seen_words.size(); ++w)
        {
          const Word wd = seen_words[w];
          if (new_word.is_equal(wd))// || new_word.is_equal_inverse(wd))
          {
            seen = true;
            break;
          }
        }

        if (!seen)
          seen_words.push_back(new_word);
      }

    }
    i += 1;
  }

  return seen_words;
}

bool CheckWords(const std::vector<Word>& Words,
                bool& SeenParabolic)
{
  SeenParabolic = false;
  for (size_t i = 0; i< Words.size(); ++i)
  {
    Word word_test = Words[i];

    if (word_test.is_non_finite_elliptic())
      return false;

    if (!SeenParabolic)
      if ((word_test.get_isom_class() == IsomClass::ParabolicPure) ||
          (word_test.get_isom_class() == IsomClass::ParabolicScrew))
        SeenParabolic = true;
  }
  return true;
}

double summary(const std::vector<Word>& Words)
{
  unsigned int ID_count = 0;
  unsigned int EL_count = 0; std::map<int, int> EL_ords;
  unsigned int RP_count = 0; std::map<int, int> RP_ords;
  unsigned int RL_count = 0; std::map<int, int> RL_ords;
  unsigned int PP_count = 0;
  unsigned int PS_count = 0;
  unsigned int LX_count = 0;

  for (size_t i = 0; i < Words.size(); ++i)
  {
    const Word w = Words[i];
    switch (w.get_isom_class())
    {
      case IsomClass::Identity:
      {
        ++ID_count;
        break;
      }
      case IsomClass::Elliptic:
      {
        ++EL_count;
        ++EL_ords[w.get_order()];
        break;
      }
      case IsomClass::ReflectionPoint:
      {
        ++RP_count;
        ++RP_ords[w.get_order()];
        break;
      }
      case IsomClass::ReflectionLine:
      {
        ++RL_count;
        ++RL_ords[w.get_order()];
        break;
      }
      case IsomClass::ParabolicPure:
      {
        ++PP_count;
        break;
      }
      case IsomClass::ParabolicScrew:
      {
        ++PS_count;
        break;
      }
      case IsomClass::Loxodromic:
      {
        ++LX_count;
        break;
      }
    }
  }
  //

  std::map<int, int>::iterator it;

  std::cout << "Total words:         " << Words.size() << std::endl;
  std::cout << "Identity:            " << ID_count << std::endl;
  std::cout << "Elliptic:            " << EL_count << std::endl;
  for ( it = EL_ords.begin(); it != EL_ords.end(); ++it )
  {
    std::cout << "\t" << it->first << ": "
              << it->second  << std::endl ;
  }
  std::cout << "Reflection in Point: " << RP_count << std::endl;
  for ( it = RP_ords.begin(); it != RP_ords.end(); ++it )
  {
    std::cout << "\t" << it->first << ": "
              << it->second  << std::endl ;
  }
  std::cout << "Reflection in Line:  " << RL_count << std::endl;
  for ( it = RL_ords.begin(); it != RL_ords.end(); ++it )
  {
    std::cout << "\t" << it->first << ": "
              << it->second  << std::endl ;
  }
  std::cout << "Pure Parabolic:      " << PP_count << std::endl;
  std::cout << "Screw Parabolic:     " << PS_count << std::endl;
  std::cout << "Loxodromic:          " << LX_count << std::endl;

  std::cout << "Non-Loxodromic ratio:"
            << 1 - ((1.0 * LX_count) / (1.0 * Words.size()))
            << std::endl;

  return (1.0 * LX_count) / (1.0 * Words.size());
}