#pragma once
#include <string>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <memory>

namespace Util {

	struct Position {
		int _startRef;
		int _endRef;
		int _startTarget;
		int _endTarget;

		~Position() {};
	};

	struct Kmer {
		std::string _kmer;
		int _position;

		~Kmer() {};
	};
	//Key -> hashcode, value -> kmer and its position in segment
	std::unordered_map<int, std::vector< std::unique_ptr<Kmer>> > H;

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

} //end Util.h