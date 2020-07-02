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

comp_d parseComplex(const std::string cString)
{
  //real[+/-]Iimag
  int from = 0;
  int pm = cString.find_first_of ("+-", from + 1);
  char c = cString.at(pm);
  const double r1 = std::stod(cString.substr(from, pm - from));
  from = pm + 2;
  const double i1 = std::stod(cString.substr(from, 100));
  const comp_d m(r1, (c == '+') ? i1 : (-1 * i1));

  return m;
}

std::tuple<comp_d, comp_d, comp_d> parseMatrixRow(const std::string row)
{
  int from = 0;
  int comma = row.find (",", from);
  const comp_d m1 = parseComplex(row.substr(from, comma - from));

  from = comma + 1;
  comma = row.find (",", from);
  const comp_d m2 = parseComplex(row.substr(from, comma - from));

  from = comma + 1;
  const comp_d m3 = parseComplex(row.substr(from, 100));

  std::tuple<comp_d, comp_d, comp_d> three(m1, m2, m3);

  return three;
}

CompMat3 parseMatrix(std::ifstream& dataStream, std::string& name)
{
  std::string line;
  getline(dataStream, line);
  // line should be of the form STRING = \n
  const int equals = line.find_first_of (" =", 0);
  name = line.substr(0, equals);
  // we don't do anything with this line

  getline(dataStream, line);
  // line should be of the form "[[123+123,123+123,123+123],""
  auto [ i11, i12, i13 ] = parseMatrixRow(line.substr(2, line.size() - 4));

  getline(dataStream, line);
  // line should be of the form " [123+123,123+123,123+123],""
  auto [ i21, i22, i23 ] = parseMatrixRow(line.substr(2, line.size() - 4));

  getline(dataStream, line);
  // line should be of the form " [123+123,123+123,123+123]]"
  auto [ i31, i32, i33 ] = parseMatrixRow(line.substr(2, line.size() - 4));

  CompMat3 mat(3, 3);
  mat << i11 << i12 << i13 << arma::endr
      << i21 << i22 << i23 << arma::endr
      << i31 << i32 << i33 << arma::endr;

  return mat;
}

/*
H =
[[123+i123,123+i123,123+i123],
 [123+i123,123+i123,123+i123],
 [123+i123,123+i123,123+i123]]

R1 =
[[123+I123,123+I123,123+I123],
 [123+I123,123+I123,123+I123],
 [123+I123,123+I123,123+I123]]

R2 =
[[123+I123,123+I123,123+I123],
 [123+I123,123+I123,123+I123],
 [123+I123,123+I123,123+I123]]

R3 =
[[123+I123,123+I123,123+I123],
 [123+I123,123+I123,123+I123],
 [123+I123,123+I123,123+I123]]
*/

// returns (H, [R1, R2, R3])
std::tuple<CompMat3, std::vector<CompMat3>>  parseMatrixFile(const std::string filename,
                                                             const bool verbose)
{
  std::ifstream data_file(filename);
  std::string name = "";

  const CompMat3 H = parseMatrix(data_file, name);

  std::vector<CompMat3> gens;

  std::string line; // not used
  while (getline(data_file, line))
  {
    CompMat3 temp = parseMatrix(data_file, name);
    // check that these matrices are isometries wrt the Hermitian form
    const CompMat3 iso_check = temp.t() * H * temp;

    if (!arma::approx_equal(H, iso_check, "absdiff", TOL))
      std::cout << "Warning: Matrix " << name << " is a not an isometry." << std::endl;

    gens.push_back(temp);
  }

  std::tuple<CompMat3, std::vector<CompMat3>> group(H, gens);

  return group;
}