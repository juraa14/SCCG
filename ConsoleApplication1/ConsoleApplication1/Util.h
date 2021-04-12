#pragma once
#include <string>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <memory>
#include <functional>

namespace Util {

	struct Position {
		size_t _startRef;
		size_t _endRef;
		size_t _startTarget;
		size_t _endTarget;

		~Position() {};
	};

	struct Kmer {
		std::string _kmer;
		size_t _position;

		Kmer(std::string _kmer, int _position): _kmer(_kmer), _position(_position) {}

		~Kmer() {};
	};

	//Key -> hashcode, value -> [kmer and its position in segment]
	std::unordered_map<size_t, std::vector< std::unique_ptr<Kmer>> > local_H;

	std::string readFASTA(std::string& fileName) {

		std::ifstream inputFile(fileName);

		if (!inputFile.good()) {
			throw std::invalid_argument("Wrong file path");
		}

		std::string data, line;

		std::getline(inputFile, line);

		while (std::getline(inputFile, line)) {
			data += line;
		}

		return data;
	}

	std::vector<std::unique_ptr<Position>> recordLowercasePositions(std::string& sequence) {
		std::vector<std::unique_ptr<Position>> positions;

		int start = 0, end = 0;
		bool first_lower = true;

		for (size_t i = 0; i < sequence.size(); i++) {
			if (islower(sequence[i])) {
				if (first_lower) {
					start = i;
					first_lower = false;
				} 
				end += 1;
			}
			else {
				if (!first_lower) {
					auto pos = std::make_unique<Position>();
					pos->_startTarget = start;
					pos->_endTarget = end - 1;
					positions.push_back(std::move(pos));
				}
				start = 0;
				end = i + 1;
				first_lower = true;
			}
		}

		if (!first_lower) {
			auto pos = std::make_unique<Position>();
			pos->_startTarget = start;
			pos->_endTarget = end - 1;
			positions.push_back(std::move(pos));
		}

		return positions;
	} 

	void writePositionsToFile(std::string fileName, std::vector< std::unique_ptr<Position> >& positions) {
		std::ofstream outputFile(fileName);

		if (!outputFile.good()) {
			throw std::invalid_argument("Can't open file!");
		}

		for (auto const& pos : positions) {
			outputFile << pos->_startTarget << " " << pos->_endTarget << " " << std::endl;
		}

		outputFile.close();
	}

	void buildLocalHashTable(std::string refSegment, int kmerLength) {
		int L = refSegment.size();

		for (size_t i = 0; i < L - kmerLength + 1; i++) {
			std::string kmerStr = refSegment.substr(i, kmerLength);

			auto kmer = std::make_unique<Kmer>(kmerStr, i);

			size_t hash = std::hash<std::string>{}(kmerStr);

			local_H[hash].push_back(std::move(kmer));			
		}
	}

} //end Util.h