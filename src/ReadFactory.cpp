/* 
 * File:   ReadFactory.cpp
 * Author: veraalva
 * 
 * Created on May 6, 2016, 10:49 AM
 */
#include <dirent.h>

#include <iostream>
#include <fstream>
#include <memory>
#include <ctime>
#include <set>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>

#include "api/BamReader.h"

#include "Global.h"
#include "Exceptions.h"
#include "TimeUtils.h"
#include "bstring.h"
#include "TextParser.h"
#include "GenomeFactory.h"
#include "ReadFactory.h"

using namespace std;
using namespace parsers;
using namespace ngs;
using namespace BamTools;
using namespace genome;

ReadFactory::ReadFactory() {
}

ReadFactory::~ReadFactory() {
}

void ReadFactory::processReadAtGenomeLevel(std::string chrName, std::string sampleName, unsigned int start, unsigned int end) {
    bool done = false;
    GeneMultiSetNGS::iterator geneIt;
    try {
        geneIt = genomeFactory.findGeneUpperBound(chrName, start, end);
        for (auto it = geneIt;; --it) {
            this->processReadAtGeneLevel(*it, sampleName, start, end);
            if (done) break;
            if ((*it)->getEnd() < start) done = true;
            if (it == genomeFactory.getCurrentChr()->getGenes().begin()) break;
        }
    } catch (exceptions::NotFoundException) {
    }
}

void ReadFactory::processReadAtGeneLevel(SPtrGeneNGS gene, std::string sampleName, unsigned int start, unsigned int end) {
    if (!gene->isInside(start, end, 8)) return;
    if (!gene->isProcessed()) gene->setProcessed(true);
    for (auto it = gene->getIsoforms().begin(); it != gene->getIsoforms().end(); ++it) {
        this->processReadAtIsoformLevel(*it, sampleName, start, end);
    }
}

void ReadFactory::processReadAtIsoformLevel(SPtrIsoformNGS isoform, std::string sampleName, unsigned int start, unsigned int end) {
    if (!isoform->isInside(start, end, 8)) return;
    isoform->setProcessed(true);
    SPtrSampleData s = isoform->getData().createSampleData(sampleName);
    s->increaseReads();
    for (auto it = isoform->getFeatures().begin(); it != isoform->getFeatures().end(); ++it) {
        if ((*it)->isInside(start, end, 8)) {
            (*it)->getData().increaseReads(sampleName);
            if (end > (*it)->getEnd() && it != --(isoform->getFeatures().end())) {
                s->increaseBridgeReads();
            }
        }
    }
}

void ReadFactory::PopulateReads(std::string sampleName) {
    SPtrChromosomeNGS c;
    SPtrGeneNGS g;
    SPtrIsoformNGS i;
    SPtrFeatureNGS f;

    for (auto cIt = genomeFactory.getChromosomes().begin(); cIt != genomeFactory.getChromosomes().end(); ++cIt) {
        c = cIt->second;
        for (auto gIt = c->getGenes().begin(); gIt != c->getGenes().end(); ++gIt) {
            g = *gIt;
            if (g->isProcessed()) {
                for (auto iIt = g->getIsoforms().begin(); iIt != g->getIsoforms().end(); ++iIt) {
                    i = *iIt;
                    SPtrSampleData s = i->getData().createSampleData(sampleName);
                    for (auto fIt = i->getFeatures().begin(); fIt != i->getFeatures().end(); ++fIt) {
                        f = *fIt;
                        try {
                            if (f->getData().getSampleData(sampleName)->getReads() != 0 && f->getType().compare("exon") == 0) {
                                s->increaseReads(f->getData().getSampleData(sampleName)->getReads());
                            }
                        } catch (exceptions::NotFoundException) {
                        }
                    }
                }
            }
        }
    }
}

void ReadFactory::calculateTPMperSample(std::string sampleName) {
    SPtrChromosomeNGS c;
    SPtrGeneNGS g;
    SPtrIsoformNGS i;
    SPtrFeatureNGS f;
    double T = 0.0;
    double TIsoformExon = 0.0;
    double TIsoformIntron = 0.0;
    double TBridges = 0.0;


    for (auto cIt = genomeFactory.getChromosomes().begin(); cIt != genomeFactory.getChromosomes().end(); ++cIt) {
        c = cIt->second;
        for (auto gIt = c->getGenes().begin(); gIt != c->getGenes().end(); ++gIt) {
            g = *gIt;
            if (g->isProcessed()) {
                for (auto iIt = g->getIsoforms().begin(); iIt != g->getIsoforms().end(); ++iIt) {
                    i = *iIt;
                    if (i->isProcessed()) {

                        try {
                            SPtrSampleData s = i->getData().createSampleData(sampleName);

                            T += static_cast<double> (s->getReads()) / static_cast<double> (i->getLength());
                            if (i->getFeatures().size() > 1) {
                                TBridges += static_cast<double> (s->getBridgeReads()) / static_cast<double> ((i->getFeatures().size() - 1));
                            }
                            for (auto fIt = i->getFeatures().begin(); fIt != i->getFeatures().end(); ++fIt) {
                                f = *fIt;
                                try {
                                    SPtrSampleData sf = f->getData().createSampleData(sampleName);

                                    if (f->getType().compare("exon") == 0) {
                                        s->increaseExonReads(sf->getReads());
                                        s->increaseExonLength(f->getLength());
                                    } else if (f->getType().compare("intron") == 0) {
                                        s->increaseIntronReads(sf->getReads());
                                        s->increaseIntronLength(f->getLength());
                                    }
                                } catch (exceptions::NotFoundException) {
                                }
                            }
                            if (s->getExonLength() > 0) {
                                TIsoformExon += static_cast<double> (s->getExonReads()) / static_cast<double> (s->getExonLength());
                            }
                            if (s->getIntronLength() > 0) {
                                TIsoformIntron += static_cast<double> (s->getIntronReads()) / static_cast<double> (s->getIntronLength());
                            }
                        } catch (exceptions::NotFoundException) {
                        }
                    }
                }
            }
        }
    }

    for (auto cIt = genomeFactory.getChromosomes().begin(); cIt != genomeFactory.getChromosomes().end(); ++cIt) {
        c = cIt->second;
        for (auto gIt = c->getGenes().begin(); gIt != c->getGenes().end(); ++gIt) {
            g = *gIt;
            if (g->isProcessed()) {
                for (auto iIt = g->getIsoforms().begin(); iIt != g->getIsoforms().end(); ++iIt) {
                    i = *iIt;
                    if (i->isProcessed()) {
                        try {
                            SPtrSampleData s = i->getData().getSampleData(sampleName);
                            s->setTPM(static_cast<double> (s->getReads() * 1.0E6) / static_cast<double> (i->getLength() * T));

                            if (i->getFeatures().size() > 1) {
                                s->setTPMBridges(static_cast<double> (s->getBridgeReads() * 1.0E6) / static_cast<double> ((i->getFeatures().size() - 1) * TBridges));
                            }

                            if (s->getExonLength() != 0) {
                                s->setTPMExon(static_cast<double> (s->getExonReads() * 1.0E6) / static_cast<double> (s->getExonLength() * TIsoformExon));
                            }
                            if (s->getIntronLength() != 0) {
                                s->setTPMIntron(static_cast<double> (s->getIntronReads() * 1.0E6) / static_cast<double> (s->getIntronLength() * TIsoformIntron));
                            }
                            for (auto fIt = i->getFeatures().begin(); fIt != i->getFeatures().end(); ++fIt) {
                                f = *fIt;
                                try {
                                    SPtrSampleData sf = f->getData().createSampleData(sampleName);

                                    if (f->getType().compare("exon") == 0) {
                                        sf->setTPM(static_cast<double> (sf->getReads() * 1.0E6) / static_cast<double> (f->getLength() * TIsoformExon));
                                    } else if (f->getType().compare("intron") == 0) {
                                        sf->setTPM(static_cast<double> (sf->getReads() * 1.0E6) / static_cast<double> (f->getLength() * TIsoformIntron));
                                    }
                                } catch (exceptions::NotFoundException) {
                                }
                            }
                        } catch (exceptions::NotFoundException(ex)) {
                        }
                    }
                }

            }
        }
    }
}

int ReadFactory::processReadsFromBAM(std::string bamFileName, std::string sampleName) {
    BamAlignment al;
    CigarOp c;
    bool toRun;
    int32_t start, end, len;
    int count = 0;
    if (!reader.Open(bamFileName)) {
        throw exceptions::NotFoundException("Can't open BAM file with name " + bamFileName);
    }
    header = reader.GetHeader();
    references = reader.GetReferenceData();

    while (reader.GetNextAlignment(al)) {
        if (al.IsMapped() && al.IsProperPair()) {
            start = end = al.Position;
            toRun = false;
            for (auto it = al.CigarData.begin(); it != al.CigarData.end(); ++it) {
                c = *it;
                if (c.Type == 'N') {
                    toRun = true;
                    break;
                }
            }
            if (toRun) {
                len = 0;
                start = al.Position;
                for (auto it = al.CigarData.begin(); it != al.CigarData.end(); ++it) {
                    c = *it;
                    len += c.Length;
                    if (c.Type == 'N') {
                        if (len >= 8) {
                            processReadAtGenomeLevel(references[al.RefID].RefName, sampleName, start, end);
                        }
                        start = end + 1 + c.Length;
                        len = 0;
                    } else {
                        end = start + len - 1;
                    }
                }
                if (len >= 8) {
                    processReadAtGenomeLevel(references[al.RefID].RefName, sampleName, start, end);
                }
            } else {
                processReadAtGenomeLevel(references[al.RefID].RefName, sampleName, al.Position, al.GetEndPosition(true, true));
            }
            count++;
        }
    }

    reader.Close();

//    PopulateReads(sampleName);
    calculateTPMperSample(sampleName);
    return count;
}

int ReadFactory::processReadsFromSAM(std::string samFileName, std::string sampleName) {
    TextParser fParser;
    int count = 0;
    CigarOp c;
    bool toRun;
    int32_t start, end, len;
    std::vector<CigarOp> cigarVec;

    try {
        fParser.setFileToParse((samFileName));
        while (fParser.iterate("#", "\t")) {
            if (fParser.getWords().size() < 9) {
                cerr << "SAM file with wrong number of fields. It should be 12 tab separated fields" << endl;
                exit(-1);
            }
            cigarVec = processCigar(fParser.getWords()[5]);

            start = end = atoi(fParser.getWords()[3].c_str());
            toRun = false;
            for (auto it = cigarVec.begin(); it != cigarVec.end(); ++it) {
                c = *it;
                if (c.Type == 'N') {
                    toRun = true;
                    break;
                }
            }
            if (toRun) {
                len = 0;
                for (auto it = cigarVec.begin(); it != cigarVec.end(); ++it) {
                    c = *it;
                    len += c.Length;
                    if (c.Type == 'N') {
                        processReadAtGenomeLevel(fParser.getWords()[2], sampleName, start, end);
                        start = end + 1 + c.Length;
                        len = 0;
                    } else {
                        end = start + len - 1;
                    }
                }
                if (len != 0) {
                    processReadAtGenomeLevel(fParser.getWords()[2], sampleName, start, end);
                }
            } else {
                processReadAtGenomeLevel(fParser.getWords()[2], sampleName, atoi(fParser.getWords()[3].c_str()), atoi(fParser.getWords()[3].c_str()) + strlen(fParser.getWords()[9].c_str()) - 1);
            }
            count++;
        }
    } catch (exceptions::FileHandledException ex) {
        std::cerr << ex.what() << std::endl;
        std::cerr << "Error parsing file" << std::endl;
        exit(-1);
    }

    PopulateReads(sampleName);
    calculateTPMperSample(sampleName);
    return count;
}

std::vector<CigarOp> ReadFactory::processCigar(std::string cigar) {
    vector<CigarOp> cigarVec;
    string value;

    for (auto it = cigar.begin(); it != cigar.end(); ++it) {
        if (std::isdigit(*it)) {
            value += (*it);
        } else {
            CigarOp v((*it), atoi(value.c_str()));
            cigarVec.push_back(v);
            value.erase();
        }
    }
    return cigarVec;
}

int ReadFactory::processBAMSAMFromDir(std::string dirName, std::vector<std::string> groupFeatures) {
    int count = 0;
    int totalCount = 0;
    struct dirent *dp;
    string sufix;
    string BAMsufix(".bam");
    string SAMsufix(".sam");
    TimeUtils uTime;
    string sampleFileName;

    if (groupFeatures.size() != 2) {
        cerr << "Sorry, this program can only manage 2 groups of data" << endl;
    }
    DIR *dirp = (DIR *) opendir(dirName.c_str());
    if (!dirp) {
        cerr << "Can't open directory: " << dirName << endl;
        exit(-1);
    }

    while ((dp = readdir(dirp)) != NULL) {
        sufix.clear();
        string fName(dp->d_name);
        if (fName[0] != '.' && fName.size() >= 4) {
            if (fName.compare(fName.size() - BAMsufix.size(), BAMsufix.size(), BAMsufix) == 0) sufix = BAMsufix;
            else if (fName.compare(fName.size() - BAMsufix.size(), BAMsufix.size(), BAMsufix) == 0)sufix = SAMsufix;
        }
        if (!sufix.empty()) {
            string fileName = dirName + "/" + fName;
            sampleFileName = fName;
            sampleFileName.replace(fName.size() - sufix.size(), sufix.size(), "");
            uTime.setTime();
            cerr << "Processing sample: " << sampleFileName;
            samples.push_back(sampleFileName);
            cerr.flush();
            try {
                if (sufix.compare(BAMsufix) == 0) {
                    count = processReadsFromBAM(fileName, sampleFileName);
                } else {
                    count = processReadsFromSAM(fileName, sampleFileName);
                }
                totalCount += count;
            } catch (exceptions::NotFoundException) {
                cerr << "Can't open file " << fileName << endl;
                exit(-1);
            }
            cerr << " in " << uTime.getElapseTimeSec() << " seconds. Processed " << count << " reads" << endl;
        }
    }
    closedir(dirp);

    return totalCount;
}

void ReadFactory::printResults(std::vector<std::string> groupFeatures) {
    int count = 0;
    string sampleFileName;
    SPtrChromosomeNGS c;
    SPtrGeneNGS g;
    SPtrIsoformNGS i;
    SPtrFeatureNGS f;
    ofstream out, ent, exp, tpm;

    for (auto sIt = samples.begin(); sIt != samples.end(); ++sIt) {
        sampleFileName = *sIt + ".out";
        out.open(sampleFileName);
        if (!out.is_open()) {
            cerr << "Can't open file " << sampleFileName << endl;
            exit(-1);
        }
        out << "Gene_Id\tTranscript_Id\tChr\tLength\tCount_Reads\tTPM\tExon_Length\tExon_Count_Reads\tExon_TPM\tIntron_Length\tIntron_Count_Reads\tIntron_TPM" << endl;

        sampleFileName = *sIt + ".ent";
        ent.open(sampleFileName);
        ent << "Gene_Id\tTranscript_Id\tChr\tType\tType_Number\tstart\tend\tLength\tCount_Reads\tTPM" << endl;
        if (!ent.is_open()) {
            cerr << "Can't open file " << sampleFileName << endl;
            exit(-1);
        }

        for (auto cIt = genomeFactory.getChromosomes().begin(); cIt != genomeFactory.getChromosomes().end(); ++cIt) {
            c = cIt->second;
            for (auto it = c->getGenes().begin(); it != c->getGenes().end(); ++it) {
                g = *it;
                if (g->isProcessed()) {
                    for (auto it1 = g->getIsoforms().begin(); it1 != g->getIsoforms().end(); ++it1) {
                        i = *it1;
                        if (i->isProcessed()) {
                            out << g->getId()
                                    << "\t" << i->getId()
                                    << "\t" << c->getId()
                                    << "\t" << i->getLength();
                            try {
                                SPtrSampleData s = i->getData().getSampleData(*sIt);                                
                                out << "\t" << s->getReads()
                                        << "\t" << s->getTPM()
                                        << "\t" << s->getExonLength()
                                        << "\t" << s->getExonReads()
                                        << "\t" << s->getTPMExon()
                                        << "\t" << s->getIntronLength()
                                        << "\t" << s->getIntronReads()
                                        << "\t" << s->getTPMIntron()
                                        << endl;
                            } catch (exceptions::NotFoundException) {
                                out << "\t0\t0.000000\t0\t0\t0.000000\t0\t0\t0.000000" << endl;
                            }
                            count = 1;
                            for (auto fIt = i->getFeatures().begin(); fIt != i->getFeatures().end(); ++fIt) {
                                f = *fIt;
                                try {
                                    SPtrSampleData s = f->getData().getSampleData(*sIt);
                                    ent << g->getId()
                                            << "\t" << i->getId()
                                            << "\t" << c->getId()
                                            << "\t" << f->getType()
                                            << "\t" << count++
                                            << "\t" << f->getStart()
                                            << "\t" << f->getEnd()
                                            << "\t" << f->getLength()
                                            << "\t" << s->getReads()
                                            << "\t" << s->getTPM()
                                            << endl;
                                } catch (exceptions::NotFoundException) {
                                    ent << g->getId()
                                            << "\t" << i->getId()
                                            << "\t" << c->getId()
                                            << "\t" << f->getType()
                                            << "\t" << count++
                                            << "\t" << f->getStart()
                                            << "\t" << f->getEnd()
                                            << "\t" << f->getLength()
                                            << "\t0\t0.000000"
                                            << endl;
                                }
                            }
                        }
                    }
                }
            }
        }
        out.close();
        ent.close();
    }

    cout << "Printing mean values" << endl;
    sampleFileName = "intron_count_per_samples.txt";
    out.open(sampleFileName);
    if (!out.is_open()) {
        cerr << "Can't open file Intron_count_per_samples.txt" << endl;
        exit(-1);
    }

    sampleFileName = "transcript_count_per_samples.txt";
    exp.open(sampleFileName);
    if (!exp.is_open()) {
        cerr << "Can't open file Transcript_expresion.txt" << endl;
        exit(-1);
    }

    sampleFileName = "transcript_all_per_sample.txt";
    tpm.open(sampleFileName);
    if (!tpm.is_open()) {
        cerr << "Can't open file transcript_all_per_sample.txt" << endl;
        exit(-1);
    }

    out << "Gene_Id\tTranscript_Id\tIntron_Number\tIntron_Id";
    exp << "Gene_Id\tTranscript_Id";
    tpm << "Gene_Id\tTranscript_Id";
    for (auto sIt = samples.begin(); sIt != samples.end(); ++sIt) {
        out << "\t" << *sIt;
        exp << "\t" << *sIt;
        tpm << "\t" << *sIt << "_exon\t" << *sIt << "_intron";
    }
    out << endl;
    exp << endl;
    tpm << endl;

    for (auto cIt = genomeFactory.getChromosomes().begin(); cIt != genomeFactory.getChromosomes().end(); ++cIt) {
        c = cIt->second;
        for (auto it = c->getGenes().begin(); it != c->getGenes().end(); ++it) {
            g = *it;
            if (g->isProcessed()) {
                for (auto it1 = g->getIsoforms().begin(); it1 != g->getIsoforms().end(); ++it1) {
                    i = *it1;
                    if (i->isProcessed()) {
                        exp << g->getId() << "\t" << i->getId();
                        tpm << g->getId() << "\t" << i->getId();
                        for (auto sIt = samples.begin(); sIt != samples.end(); ++sIt) {
                            try {
                                SPtrSampleData s = i->getData().getSampleData(*sIt);
                                exp << "\t" << s->getExonReads();
                                tpm << "\t" << s->getTPMExon() << "\t" << s->getTPMIntron();
                            } catch (exceptions::NotFoundException) {
                                exp << "\t0";
                                tpm << "\t0.000000\t0.000000";
                            }
                        }
                        exp << endl;
                        tpm << endl;

                        count = 1;
                        for (auto fIt = i->getFeatures().begin(); fIt != i->getFeatures().end(); ++fIt) {
                            f = *fIt;
                            if (f->getType().compare("intron") == 0) {
                                out << g->getId() << "\t"
                                        << i->getId() << "\t"
                                        << count << "\t"
                                        << i->getId() << "_" << count;
                                for (auto sIt = samples.begin(); sIt != samples.end(); ++sIt) {
                                    try {
                                        SPtrSampleData s = f->getData().getSampleData(*sIt);
                                        out << "\t" << s->getReads();
                                    } catch (exceptions::NotFoundException) {
                                        out << "\t0";
                                    }
                                }
                                count++;
                                out << endl;
                            }
                        }
                    }
                }
            }
        }
    }
    out.close();
    exp.close();
    tpm.close();
}

int ReadFactory::processReadsFromIsoformBAM(std::string bamFileName, std::string sample) {
    BamAlignment al;
    CigarOp c;
    //    bool toRun;
    int32_t rstart, rend, rlength, length;
    int32_t start;
    //    int32_t end, rlength;
    int32_t fStart, fEnd;
    int count = 0;
    if (!reader.Open(bamFileName)) {
        throw exceptions::NotFoundException("Can't open BAM file with name " + bamFileName);
    }
    header = reader.GetHeader();
    references = reader.GetReferenceData();

    while (reader.GetNextAlignment(al)) {
        bool p = references.at(al.RefID).RefName.compare("NR_027232") == 0 ? true : false;
        p = false;

        if (p) {
            cout << references.at(al.RefID).RefName << "\t" << count << endl;
        }
        try {
            SPtrIsoformNGS i = genomeFactory.findIsoform(references.at(al.RefID).RefName);
            SPtrGeneNGS g = genomeFactory.getCurrentGene();

            rlength = al.Length;
            rstart = al.Position;
            rend = al.GetEndPosition(true, false);

            fStart = 0;
            fEnd = 0;

            g->setProcessed(true);
            i->setProcessed(true);
            i->getData().increaseReads(sample);
            SPtrSampleData s = i->getData().getSampleData(sample);
            if (p) {
                cout << i->getId() << "\t" << s->getReads() << "\t" << s->getBridgeReads() << endl;
            }
            for (auto fIt = i->getFeatures().begin(); fIt != i->getFeatures().end(); ++fIt) {
                SPtrFeatureNGS f = *fIt;

                start = f->getStart();
                length = f->getEnd() - f->getStart() + 1;

                //                if (i->getStrand() == '-' && fIt == --(i->getFeatures().end())) {
                //                    length -= 1;
                //                }

                fEnd = fStart + length - 1;
                if (p) {
                    cout << "\t" << rstart << "\t" << rend << "\t" << fStart << "\t" << fEnd << "\t" << length;
                }
                if (fIt != i->getFeatures().begin() &&
                        rstart <= fStart && rend >= fStart) {
                    s->increaseBridgeReads();
                }
                if (p) {
                    cout << "\t" << s->getBridgeReads() << endl;
                }
                fStart = fEnd + 1;
            }
            count++;
        } catch (exceptions::NotFoundException) {
            if (p) {
                cout << "Error" << endl;
            }
        }
    }
    calculateTPMperSample(sample);
    reader.Close();
    return count;
}

