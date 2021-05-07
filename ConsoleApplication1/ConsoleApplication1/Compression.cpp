#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>

#include "Util.h"


int g_k = 21;        //kmer length
int g_L = 30;  //Segment size
int g_m = 100;       //LIMIT search length when global matching
double g_T1 = 0.5;   //Threshold 1
double g_T2 = 4;     //Threshold 2
bool g_global = false; //Perform global or local matching
int g_mismatch;
int g_limit = 100; //global search kmer extension limit


std::string tempFile = "tempFile";
std::string compressedFile = "compressedFile";


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
    std::string tarSequence = Util::readFASTA(targetGenome);
    std::string refSequence = Util::readFASTA(referenceGenome);

    bool lowercase_exist = true;

    std::string refSegment, tarSegment;

    while(g_local()) {

        //Get lowercase positions and write them to file
        std::vector<std::unique_ptr<Util::Position>> Lposition_target = Util::recordLowercasePositions(tarSequence);
        int targetLength = tarSequence.size();

        //std::vector<std::unique_ptr<Util::Position>> Lposition_reference = Util::recordLowercasePositions(referentialSequence);

        g_mismatch = 0;

        if (Lposition_target.size() == 0) {
            lowercase_exist = false;
        }
        else {
            Util::writePositionsToFile("target-lowercase-position.txt", Lposition_target);
            std::transform(tarSequence.begin(), tarSequence.end(), tarSequence.begin(), std::toupper);
        }
        std::transform(refSequence.begin(), refSequence.end(), refSequence.begin(), std::toupper);

        std::ofstream outTempFile(tempFile, std::ios::app);

        //kmer search
        int startRef = 0;
        int endRef = g_L;
        int startTar = 0;
        int endTar = g_L;
        auto pos = std::make_unique<Util::Position>();

        while (true) {
            
            //Util::buildLocalHashTable("ABCCAB", 2);

            if (endRef > refSequence.size()) {
                std::string restOfRefSequence = tarSequence.substr(startRef);
                outTempFile << restOfRefSequence;
                break;
            }
            if (endTar > tarSequence.size()) {
                std::string restOfRefSequence = tarSequence.substr(startRef);
                outTempFile << restOfRefSequence;
                break;
            }

            refSegment = refSequence.substr(startRef, endRef - startRef);
            tarSegment = tarSequence.substr(startTar, endTar - startTar);

            auto positions = Util::localMatching("ABADABAcABA", "ABACABAxABA", 2);
            //auto position2 = Util::localMatching("CGGACGC", "CGGGGACGGCA", 2);
            if (positions.size() == 0) {
                g_k = 11;
                positions = Util::localMatching("ABADABAA", "ABACABAA", 2);
            }

            if (positions.size() == 0) { // No local match
                g_mismatch++;

                outTempFile << tarSegment << "\n";

                startTar += g_L;
                endRef += g_L;
                endTar = startRef + g_L;

                int tar_charsTillEnd = tarSegment.size() - startTar;
                int ref_charsTillEnd = refSegment.size() - startRef;

                if (tar_charsTillEnd <= g_k) { // num of chars left in target to process is less than kmer size
                    outTempFile << tarSequence.substr(startTar);
                    break;
                }
                else if (tar_charsTillEnd < g_L) { // num of chars left in target to process is less than segment size
                    endTar = tarSequence.size() - 1;
                }

                if (ref_charsTillEnd < g_L) {
                    endRef = refSequence.size() - 1;
                }

                if (endTar >= tarSequence.size()) {
                    break;
                }

                if (ref_charsTillEnd <= g_k) {
                    std::string temp = tarSequence.substr(startTar);
                    if (temp.empty()) {
                        break;
                    }
                    else {
                        if (temp.size() > g_L * g_T1) {
                            g_mismatch++;
                        }

                        if (g_mismatch > g_T2) {
                            g_global = true;
                            break;
                        }
                        outTempFile << temp;
                        break;
                    }
                }

                continue;
            }

            break;
        }

        break;
    }

    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "Compression time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;

    return 0;
}

