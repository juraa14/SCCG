#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>

#include "Util.h"


int g_k = 21;        //kmer length
int g_L = 30;  //Segment size
int g_m = 100;       //LIMIT
double g_T1 = 0.5;   //Threshold 1
double g_T2 = 4;     //Treshold 2
bool g_global = false; //Perform global or local matching
int g_mismatch;

inline bool g_local() {
    return !g_global;
}


int main(int argc, char **argv)
{   
    if (argc != 3) {
        std::cout << "Wrong number of arguments!" << std::endl;
        std::exit(-1);
    }
    auto start = std::chrono::high_resolution_clock::now();

    std::string referenceGenome = std::string(argv[2]);
    std::string targetGenome = std::string(argv[1]);

    std::cout << "Compressing " + targetGenome << std::endl;

    //Read sequences
    std::string targetSequence = Util::readFASTA(targetGenome);
    std::string referentialSequence = Util::readFASTA(referenceGenome);

    bool lowercase_exist = true;

    std::string refSegment, targetSegment;

    while(g_local()) {

        //Get lowercase positions and write them to file
        std::vector<std::unique_ptr<Util::Position>> Lposition_target = Util::recordLowercasePositions(targetSequence);
        int targetLength = targetSequence.size();

        //std::vector<std::unique_ptr<Util::Position>> Lposition_reference = Util::recordLowercasePositions(referentialSequence);

        g_mismatch = 0;

        if (Lposition_target.size() == 0) {
            lowercase_exist = false;
        }
        else {
            Util::writePositionsToFile("target-lowercase-position.txt", Lposition_target);
            std::transform(targetSequence.begin(), targetSequence.end(), targetSequence.begin(), std::toupper);
        }
        std::transform(referentialSequence.begin(), referentialSequence.end(), referentialSequence.begin(), std::toupper);



        //kmer search
        int startRef = 0;
        int endRef = g_L;
        int startTarget = 0;
        int endTarget = g_L;
        auto pos = std::make_unique<Util::Position>();

        while (true) {

            refSegment = referentialSequence.substr(startRef, endRef - startRef);
            targetSegment = targetSequence.substr(startTarget, endTarget - startTarget);

            int k = g_k; //kmer length k
            


            break;
        }

        break;
    }

    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "Compression time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;

    return 0;
}

