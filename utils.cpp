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

std::string GenerateFileName()
{
  time_t theTime = time(NULL);
  struct tm *aTime = localtime(&theTime);

  int day   = aTime->tm_mday;
  int month = aTime->tm_mon + 1;
  int year  = aTime->tm_year + 1900;
  int hour  = aTime->tm_hour;
  int min   = aTime->tm_min;
  int sec   = aTime->tm_sec;

  std::string filename = "output/";
  char buffer [50];
  sprintf(buffer, "%04d", year);
  filename.append(buffer);
  sprintf(buffer, "%02d", month);
  filename.append(buffer);
  sprintf(buffer, "%02d", day);
  filename.append(buffer);
  filename.append("-");
  sprintf(buffer, "%02d", hour);
  filename.append(buffer);
  sprintf(buffer, "%02d", min);
  filename.append(buffer);
  sprintf(buffer, "%02d", sec);
  filename.append(buffer);
  filename.append(".txt");

  return filename;
}