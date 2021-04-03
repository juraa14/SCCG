#pragma once
#include <string>
#include <fstream>
#include <vector>

namespace Util {

	struct Position {
		int startRef;
		int endRef;
		int startTarget;
		int endTarget;
	};

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

	std::vector<Position*> recordLowercasePositions(std::string sequence) {
		std::vector<Position*> positions;

		int start = 0, end = 0;
		bool first_lower = true;

		for (size_t i = 0; i < sequence.size(); i++) {
			if (islower(sequence[i])) {
				if (first_lower) {
					start = i;
					end += 1;
					first_lower = false;
				} 
				else {
					end += 1;
				}
			}
			else {
				if (!first_lower) {
					Position* pos = new Position();
					pos->startTarget = start;
					pos->endTarget = end - 1;
					positions.push_back(pos);
				}
				start = 0;
				end = i + 1;
				first_lower = true;
			}
		}

		if (!first_lower) {
			Position* pos = new Position();
			pos->startTarget = start;
			pos->endTarget = end - 1;
			positions.push_back(pos);
		}

		return positions;
	}

} //end Util.h