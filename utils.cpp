#include "includes.h"
#include "utils.h"

unsigned int gcd(const unsigned int m, const unsigned int n)
{
  unsigned int r = m % n;
  if (r != 0)
    return gcd(n, r);
  else
    return n;
}
