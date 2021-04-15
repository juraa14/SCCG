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

	std::vector<std::unique_ptr<Position>> recordLowercasePositions(const std::string& sequence) {
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

	void writePositionsToFile(const std::string& fileName, std::vector< std::unique_ptr<Position> >& positions) {
		std::ofstream outputFile(fileName);

		if (!outputFile.good()) {
			throw std::invalid_argument("Can't open file!");
		}

		for (auto const& pos : positions) {
			outputFile << pos->_startTarget << " " << pos->_endTarget << " " << std::endl;
		}

		outputFile.close();
	}

	void buildLocalHashTable(const std::string& refSegment, int kmerLength) {
		size_t L = refSegment.size();

		for (size_t i = 0; i < L - kmerLength + 1; i++) {
			std::string kmerStr = refSegment.substr(i, kmerLength);

			auto kmer = std::make_unique<Kmer>(kmerStr, i);

			size_t hash = std::hash<std::string>{}(kmerStr);

			local_H[hash].push_back(std::move(kmer));			
		}
	}

	std::vector<std::unique_ptr<Position>> localMatching(const std::string& reference, const std::string& target, int kmerLength) {
		std::vector<std::unique_ptr<Position>> positions;

		int startPosition, incrementSize, longestIncrement;
		
		buildLocalHashTable(reference, kmerLength);

		for (size_t i = 0; i < reference.size() - kmerLength + 1; i++) {
			std::string targetKmer = target.substr(i, kmerLength);
			size_t hashTar = std::hash<std::string>{}(targetKmer);

			if (local_H.find(hashTar) == local_H.end()) { //target & reference segment do not match
				startPosition = -1;
				continue; //Increment index
			}

			auto& kmerVecRef = local_H[hashTar]; //Get kmer positions reference segment
			incrementSize = 0; //No kmer increments for now
			longestIncrement = 0; //If there are multiple extended segments with same length

			for (auto const& kmerRef : kmerVecRef) {
				incrementSize = 0;
				if (kmerRef->_kmer == targetKmer) {
					size_t i_ref = reference.size() - 1;
					size_t i_tar = target.size() - 1;
					size_t endRef = kmerRef->_position + kmerLength - 1;
					size_t endTar = i + kmerLength - 1;

					//How many chars can we extend the matching kmers?
					for (auto [i_endRef, i_endTar] = std::tuple{ endRef + 1, endTar + 1 }; i_endRef <= i_ref && i_endTar <= i_tar && reference[i_endRef] == target[i_endTar]; i_endRef++, i_endTar++) {
						incrementSize++;
					}

					if (kmerVecRef.size() <= 1) {
						longestIncrement = incrementSize;
						startPosition = kmerRef->_position;
						continue;
					}

					else if (incrementSize == longestIncrement) { //Found same largest increment
						if (positions.size() > 1) {
							size_t endRefLast = positions.back()->_endRef;
							if (kmerRef->_position - endRef < startPosition - endRefLast && startPosition != -1) {
								startPosition = kmerRef->_position;
							}
						}
					}
					else {
						longestIncrement = incrementSize;
						startPosition = kmerRef->_position;
					}
				}
			}

			if (startPosition == -1) { //No extended matching found
				continue;
			}

			auto pos = std::make_unique<Position>();
			pos->_startTarget = i;
			pos->_endTarget = i + kmerLength + longestIncrement - 1;
			pos->_startRef = startPosition;
			pos->_endRef = startPosition + kmerLength + longestIncrement - 1;

			positions.push_back(std::move(pos));

			i += kmerLength + longestIncrement - 2; //-2
		}

		return positions;
	}

} //end Util.h