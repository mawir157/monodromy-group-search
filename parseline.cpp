#include "parsedline.h"

ParsedLine::ParsedLine(const unsigned int a11, const unsigned int a12,
											 const unsigned int a21, const unsigned int a22,
											 const unsigned int a31, const unsigned int a32,
											 const unsigned int b11, const unsigned int b12,
											 const unsigned int b21, const unsigned int b22,
											 const unsigned int b31, const unsigned int b32,
											 const int ref_order) :
		a11(a11), a12(a12), a21(a21), a22(a22), a31(a31), a32(a32),
		b11(b11), b12(b12), b21(b21), b22(b22), b31(b31), b32(b32),
		ref_order(ref_order)
{}

std::string ParsedLine::print_helper(const unsigned int num,
			                               const unsigned int den,
			                    					 const std::string str) const
{
	std::string s;
  s.append(std::to_string(num));
  if (num != 0)
  {
	  s.append("/");
	  s.append(std::to_string(den));
	}
	s.append(str);
	return s;
}

std::string ParsedLine::print() const
{
	std::string s = "(";
	s.append(print_helper(a11, a12, ","));
	s.append(print_helper(a21, a22, ","));
	s.append(print_helper(a31, a32, ";"));
	s.append(print_helper(b11, b12, ","));
	s.append(print_helper(b21, b22, ","));
	s.append(print_helper(b31, b32, ")"));

  s.append("<");
  s.append(std::to_string(ref_order));
  s.append(">");

  return s;
}

std::tuple<CompMat3, CompMat3, CompMat3> ParsedLine::toMatrices() const
{
	const Triple ATriple(a11,a12, a21,a22, a31,a32);
	const Triple BTriple(b11,b12, b21,b22, b31,b32);

	const CompMat3 A = TripleToMatrix(ATriple);
	const CompMat3 B = TripleToMatrix(BTriple);
	const CompMat3 H = SextupleToH(ATriple, BTriple);

	std::tuple<CompMat3, CompMat3, CompMat3> three(A, B, H);

	return three;
}

std::vector<ParsedLine> parseFile(const std::string filename,
																	const bool verbose)
{
	std::string line;
	std::ifstream data_file(filename);
	std::vector<ParsedLine> pf;

	while (getline(data_file, line))
	{
		//std::cout << std::endl << line << std::endl << std::endl;
		if (line.substr(0, 1) == "#")
		{
			if (verbose)
				std::cout << line << std::endl;
			
			continue;
		}

		if (line.substr(0, 1) == "!")
		{
			continue;
		}

		if (line == "")
			continue;

		int from = 0;
		int comma = line.find(",", from);
		const int a11 = std::stoi(line.substr(from, comma - from));
		from = comma + 1;
		comma = line.find(",", from);
		const int a12 = std::stoi(line.substr(from, comma - from));
		from = comma + 1;
		comma = line.find(",", from);
		const int a21 = std::stoi(line.substr(from, comma - from));
		from = comma + 1;
		comma = line.find(",", from);
		const int a22 = std::stoi(line.substr(from, comma - from));
		from = comma + 1;
		comma = line.find(",", from);
		const int a31 = std::stoi(line.substr(from, comma - from));
		from = comma + 1;
		comma = line.find(",", from);
		const int a32 = std::stoi(line.substr(from, comma - from));

		from = comma + 1;;
		comma = line.find(",", from);
		const int b11 = std::stoi(line.substr(from, comma - from));
		from = comma + 1;
		comma = line.find(",", from);
		const int b12 = std::stoi(line.substr(from, comma - from));
		from = comma + 1;
		comma = line.find(",", from);
		const int b21 = std::stoi(line.substr(from, comma - from));
		from = comma + 1;
		comma = line.find(",", from);
		const int b22 = std::stoi(line.substr(from, comma - from));
		from = comma + 1;
		comma = line.find(",", from);
		const int b31 = std::stoi(line.substr(from, comma - from));
		from = comma + 1;
		comma = line.find(",", from);
		const int b32 = std::stoi(line.substr(from, comma - from));

		// expected reflection order
		from = comma + 1;
		comma = line.find(",", from);
		const std::string ref_order = line.substr(from, comma - from);
		int p = -1;
		if (ref_order != "*")
			p = std::stoi(ref_order);

		pf.emplace_back(a11,a12,a21,a22,a31,a32,
			              b11,b12,b21,b22,b31,b32,
			              p);
	}	
	return pf;
}

std::string GenerateFileName(const std::string directory,
														 const std::string label)
{
  time_t theTime = time(NULL);
  struct tm *aTime = localtime(&theTime);

  int day   = aTime->tm_mday;
  int month = aTime->tm_mon + 1;
  int year  = aTime->tm_year + 1900;
  int hour  = aTime->tm_hour;
  int min   = aTime->tm_min;
  int sec   = aTime->tm_sec;

  std::string filename = "./";
  filename.append(directory);
  filename.append("/");
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

  if (label.length() > 0)
		filename.append(label);

  filename.append(".txt");

  return filename;
}