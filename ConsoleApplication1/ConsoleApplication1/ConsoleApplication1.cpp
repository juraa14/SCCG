#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>

#include "Kmer.h"
#include "Util.h"


int main(int argc, char **argv)
{   
    if (argc != 3) {
        std::cout << "Wrong number of arguments!" << std::endl;
        std::exit(-1);
    }
    auto start = std::chrono::high_resolution_clock::now();

    std::string referenceGenome = std::string(argv[2]);
    std::string targetGenome = std::string(argv[1]);

    std::string targetSequence = Util::readFASTA(targetGenome);
    std::string referentialSequence = Util::readFASTA(referenceGenome);

    std::cout << referentialSequence << std::endl;

    std::string testStr = "aaaaAAaBBaa";
    std::vector<Util::Position*> pos = Util::recordLowercasePositions(testStr);

    Util::writePositionsToFile("output-Proba.txt", pos);

    std::transform(testStr.begin(), testStr.end(), testStr.begin(), std::toupper);
    std::cout << testStr << std::endl;

    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "Execution time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
}

