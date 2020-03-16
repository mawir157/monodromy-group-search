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

std::string ParsedLine::print() const
{
	std::string s = "(";
  s.append(std::to_string(a11));
  if (a11 != 0)
  {
	  s.append("/");
	  s.append(std::to_string(a12));
	}
	s.append(",");

  s.append(std::to_string(a21));
  s.append("/");
  s.append(std::to_string(a22));
  s.append(",");

  s.append(std::to_string(a31));
  s.append("/");
  s.append(std::to_string(a32));
  s.append(";");
 
 	s.append(std::to_string(b11));
  if (b11 != 0)
  {
  	s.append("/");
  	s.append(std::to_string(b12));
  }
  s.append(",");

  s.append(std::to_string(b21));
  s.append("/");
  s.append(std::to_string(b22));
  s.append(",");

  s.append(std::to_string(b31));
  s.append("/");
  s.append(std::to_string(b32));
  s.append(")");

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
