/* 
 * File:   FastaFactory.cpp
 * Author: veraalva
 * 
 * Created on February 10, 2016, 3:41 PM
 */

#include <dirent.h>
#include <inttypes.h>

#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <memory>
#include <unordered_map>
#include <utility>
#include <fstream>
#include <random>
#include <algorithm>
#include <ctime>
#include <chrono>

#include "Global.h"
#include "Exceptions.h"
#include "TimeUtils.h"
#include "TextParser.h"
#include "Sequence.h"
#include "FastaFactory.h"

using namespace std;
using namespace parsers;
using namespace sequence;
using namespace formats;

long unsigned int FastaFactory::parseDNAFastaFile(sequence::DNAContainer &seqContainer, std::string fName, bool binary) {
    int numberSeqCurrentRead = 0;
    uint64_t i, len;
    string id;
    std::shared_ptr<sequence::DNA> dna;
    pair < SPtrDNA, bool> result;

    if (binary) {
        ifstream inFile(fName, std::ifstream::binary);
        inFile.read((char *) &i, sizeof (uint64_t));
        for (unsigned long int j = 0; j < i; j++) {
            inFile.read((char *) &len, sizeof (uint64_t));
            id.resize(len);
            inFile.read(&(id[0]), len);
            result = seqContainer.addElement(id);
            if (result.second) {
                dna = result.first;
            } else {
                cerr << "Duplicated sequence ID " << id << endl;
                exit(-1);
            }
            dna->setId(id);
            inFile.read((char *) &len, sizeof (uint64_t));
            dna->getSeq().resize(len);
            inFile.read(&(dna->getSeq()[0]), len);
            numberSeqCurrentRead++;
        }
        inFile.close();
    } else {
        TextParser fParser;
        try {
            fParser.setFileToParse(fName);
            while (fParser.iterate("#")) {
                string line = fParser.getLine();
                if (fParser.lineStartWith(">")) {
                    id = line.substr(1, fParser.getLine().size() - 1);
                    result = seqContainer.addElement(id);
                    if (result.second) {
                        dna = result.first;
                    } else {
                        cerr << "Duplicated sequence ID " << id << endl;
                        exit(-1);
                    }
                    dna->setId(id);
                    numberSeqCurrentRead++;
                } else {
                    if (dna == nullptr) {
                        cerr << "Fasta file does not start with the header (>)" << endl;
                        exit(-1);
                    }
                    dna->getSeq().append(line);
                }
            }
        } catch (exceptions::FileHandledException) {
            cerr << "Error parsing file: " << fName << endl;
            exit(-1);
        } catch (ios::failure) {
            cerr << "Error parsing file: " << fName << endl;
            exit(-1);
        }
    }

    if (Global::instance()->isDebug3()) {
        cout << "\tDEBUG3 ==> " << seqContainer.size() << " sequences in the container." << endl;
        for (auto it = seqContainer.getContainer().begin(); it != seqContainer.getContainer().end(); ++it) {
            dna = it->second;
            cout << "\tDEBUG3 ==>\t\t" << dna->getId() << " with " << dna->getLength() << " bp" << endl;
        }
    }

    return numberSeqCurrentRead;
}

void FastaFactory::parseDNAFastaInDirectory(sequence::DNAContainer &seqContainer, std::string dirName, std::string prefix, std::string sufix, bool binary) {
    struct dirent *dp;
    TimeUtils tUtil;
    DIR *dirp = (DIR *) opendir(dirName.c_str());
    if (!dirp) {
        cerr << "Can't open directory: " << dirName << endl;
        exit(-1);
    }

    while ((dp = readdir(dirp)) != NULL) {
        bool read = false;
        string fName(dp->d_name);
        if (Global::instance()->isDebug3()) {
            cout << "\tDEBUG3 ==> Found file: " << fName << endl;
        }
        if (fName[0] != '.') {
            if (prefix.empty() && sufix.empty()) {
                read = true;
            } else {
                if (!prefix.empty() && sufix.empty()) {
                    if (prefix.size() <= fName.size() &&
                            fName.compare(0, prefix.size(), prefix) == 0) read = true;
                } else if (prefix.empty() && !sufix.empty()) {
                    if (sufix.size() <= fName.size() &&
                            fName.compare(fName.size() - sufix.size(), sufix.size(), sufix) == 0) read = true;
                } else if (!prefix.empty() && !sufix.empty()) {
                    if (prefix.size() <= fName.size() &&
                            sufix.size() <= fName.size() &&
                            fName.compare(0, prefix.size(), prefix) == 0 &&
                            fName.compare(fName.size() - sufix.size(), sufix.size(), sufix) == 0) read = true;
                }
            }
        }
        if (read) {
            if (Global::instance()->isInfo()) {
                tUtil.setTime();
                cout << "\tINFO ==> Parsing file: " << dirName + "/" + fName << endl;
            }
            int seqs = parseDNAFastaFile(seqContainer, dirName + "/" + fName, binary);
            if (Global::instance()->isInfo()) {
                cout << "\tINFO ==> " << seqs << " sequences read in " << tUtil.getElapseTimeSec() << " sec" << endl;
            }
        }
    }
    closedir(dirp);
}

void FastaFactory::writeDNASequencesToFile(sequence::DNAContainer &seqContainer, std::string fileName, bool binary) {
    uint64_t i, len;
    DNA *dna;

    if (Global::instance()->isDebug3()) {
        cout << "\tDEBUG3 ==> " << seqContainer.size() << " sequences in the container to write" << endl;
    }

    if (binary) {
        std::ofstream outputFile(fileName, std::ofstream::binary);
        if (!outputFile.is_open()) {
            cerr << "Can't open output file " << fileName << endl;
            exit(-1);
        }
        i = seqContainer.size();
        outputFile.write((char *) &i, sizeof (uint64_t));
        for (auto it = seqContainer.getContainer().begin(); it != seqContainer.getContainer().end(); ++it) {
            dna = it->second.get();
            if (Global::instance()->isDebug3()) {
                cout << "\tDEBUG3 ==> Writing sequence: " << dna->getId() << " with length " << dna->getLength() << endl;
            }
            len = dna->getId().size();
            outputFile.write((char *) &len, sizeof (uint64_t));
            outputFile.write(dna->getId().c_str(), len);
            len = dna->getLength();
            outputFile.write((char *) &len, sizeof (uint64_t));
            outputFile.write(dna->getSeq().c_str(), len);
        }
        outputFile.close();
    } else {
        ofstream outputFile(fileName);
        if (!outputFile.is_open()) {
            cerr << "Can't open output file " << fileName << endl;
            exit(-1);
        }
        for (auto it = seqContainer.getContainer().begin(); it != seqContainer.getContainer().end(); ++it) {
            dna = it->second.get();
            if (Global::instance()->isDebug3()) {
                cout << "\tDEBUG3 ==> Writing sequence: " << dna->getId() << " with length " << dna->getLength() << endl;
            }
            outputFile << ">" << dna->getId() << endl;
            for (i = 0; i < dna->getLength(); i += 50) {
                char t = 0;
                if (i + 50 < dna->getLength()) {
                    t = dna->getSeq()[i + 50];
                    dna->getSeq()[i + 50] = 0;
                }
                outputFile << (dna->getSeq().c_str() + i) << endl;
                if (i + 50 < dna->getLength()) {
                    dna->getSeq()[i + 50] = t;
                }
            }
        }
        outputFile.close();
    }

}
