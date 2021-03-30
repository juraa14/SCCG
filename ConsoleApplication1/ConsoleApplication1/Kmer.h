#pragma once
#include <string>

class Kmer
{
private:
	std::string kmer;
	int kmetstart;

public:
	Kmer(std::string kmer, int kmerstart);
};

