#pragma once
#include <string>

class Kmer
{
private:
	std::string kmer_;
	int kmerStart_;

public:
	Kmer(std::string kmer, int kmerStart);
};

