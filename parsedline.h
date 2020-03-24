#pragma once

#include "matrixfns.h"

class ParsedLine
{
	public:
		ParsedLine(const unsigned int a11, const unsigned int a12,
							 const unsigned int a21, const unsigned int a22,
							 const unsigned int a31, const unsigned int a32,
							 const unsigned int b11, const unsigned int b12,
							 const unsigned int b21, const unsigned int b22,
							 const unsigned int b31, const unsigned int b32,
							 const int ref_order);
		std::string print() const;
		std::tuple<CompMat3, CompMat3, CompMat3> toMatrices() const;
		int get_ref_order() const { return ref_order; }
	private:
		std::string print_helper(const unsigned int num,
														 const unsigned int den,
														 const std::string str) const;

		const unsigned int a11;
		const unsigned int a12;
		const unsigned int a21;
		const unsigned int a22;
		const unsigned int a31;
		const unsigned int a32;
		const unsigned int b11;
		const unsigned int b12;
		const unsigned int b21;
		const unsigned int b22;
		const unsigned int b31;
		const unsigned int b32;
		const int ref_order;
};

std::vector<ParsedLine> parseFile(const std::string filename,
	                                const bool verbose = false);

std::string GenerateFileName(const std::string directory = "output",
														 const std::string label = "");
