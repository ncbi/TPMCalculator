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
#include "Sequence.h"
#include "GenomeFactory.h"
#include "ReadFactory.h"
#include "RandomFactory.h"

using namespace std;
using namespace parsers;
using namespace sequence;
using namespace ngs;
using namespace BamTools;
using namespace genome;

ReadFactory::ReadFactory() {
}

ReadFactory::~ReadFactory() {
}

void ReadFactory::processReadAtGenomeLevel(std::string chrName, std::string sampleName, std::set < std::pair<unsigned int, unsigned int>, ReadFactory::coordinateLessCMP> read_coords, uint16_t minOverlap) {
    bool done = false;
    GeneMultiSetNGS::iterator geneIt;
    try {
        if (!read_coords.empty()) {
            geneIt = genomeFactory.findGeneUpperBound(chrName, read_coords.rbegin()->first, read_coords.rbegin()->second);
            for (auto it = geneIt;; --it) {
                this->processReadAtGeneLevel(*it, sampleName, read_coords, minOverlap);
                if (done) break;
                if (std::distance(it, geneIt) > 6 && (*it)->getEnd() < (*read_coords.begin()).first) done = true;
                if (it == genomeFactory.getCurrentChr()->getGenes().begin()) break;
            }
        }
    } catch (exceptions::NotFoundException) {
    }
}

void ReadFactory::processReadAtGeneLevel(SPtrGeneNGS gene, std::string sampleName, std::set < std::pair<unsigned int, unsigned int>, ReadFactory::coordinateLessCMP> read_coords, uint16_t minOverlap) {
    bool counted = false;
    bool countedExon = false;
    bool countedIntron = false;
    bool uniqueCounted = false;
    bool uniqueCountedExon = false;
    bool uniqueCountedIntron = false;
    bool countedBridge = false;

    for (auto coordIt = read_coords.begin(); coordIt != read_coords.end(); ++coordIt) {
        unsigned int start = coordIt->first;
        unsigned int end = coordIt->second;

        if (!gene->isInside(start, end, minOverlap)) return;
        if (!gene->isProcessed()) gene->setProcessed(true);
        SPtrSampleData s = gene->getData().createSampleData(sampleName);
        if (!counted) {
            s->increaseReads();
            counted = true;
        }

        // Counting a gene level for exons and introns created by the collapse of the isoforms
        for (auto it = gene->getFeatures().begin(); it != gene->getFeatures().end(); ++it) {
            SPtrFeatureNGS feature = *it;
            SPtrSampleData sF = feature->getData().createSampleData(sampleName);
            if (feature->isInside(start, end, minOverlap)) {
                sF->increaseReads();
                if (end > feature->getEnd() && it != --(gene->getFeatures().end())&& !countedBridge) {
                    s->increaseBridgesReads();
                    countedBridge = true;
                }
                if (feature->getType().compare("exon") == 0 && !countedExon) {
                    s->increaseExonReads();
                    countedExon = true;
                }
                if (feature->getType().compare("intron") == 0 && !countedIntron) {
                    s->increaseIntronReads();
                    countedIntron = true;
                }
            }
        }

        // Counting at gene level for Unique features
        for (auto it = gene->getUniqueFeatures().begin(); it != gene->getUniqueFeatures().end(); ++it) {
            SPtrFeatureNGS feature = *it;
            SPtrSampleData sF = feature->getData().createSampleData(sampleName);
            if (feature->isInside(start, end, minOverlap)) {
                if (!uniqueCounted) {
                    s->increaseUniqueReads();
                    uniqueCounted = true;
                }
                sF->increaseReads();
                if (feature->getType().compare("exon") == 0 && !uniqueCountedExon) {
                    s->increaseUniqueExonReads();
                    uniqueCountedExon = true;
                }
                if (feature->getType().compare("intron") == 0 && !uniqueCountedIntron) {
                    s->increaseUniqueIntronReads();
                    uniqueCountedIntron = true;
                }
            }
        }
    }

    // Assigning reads to isoforms and features
    for (auto it = gene->getIsoforms().begin(); it != gene->getIsoforms().end(); ++it) {
        SPtrIsoformNGS isoform = *it;
        if (!isoform->isProcessed()) isoform->setProcessed(true);
        bool countedIsoform = false;
        bool countedIsoformExon = false;
        bool countedIsoformIntron = false;
        bool countedIsoformBridge = false;
        for (auto coordIt = read_coords.begin(); coordIt != read_coords.end(); ++coordIt) {
            unsigned int start = coordIt->first;
            unsigned int end = coordIt->second;
            if (isoform->isInside(start, end, minOverlap)) {
                SPtrSampleData sI = isoform->getData().createSampleData(sampleName);
                if (!countedIsoform) {                    
                    sI->increaseReads();
                    countedIsoform = true;
                }
                for (auto it2 = isoform->getFeatures().begin(); it2 != isoform->getFeatures().end(); ++it2) {
                    SPtrFeatureNGS feature = *it2;
                    SPtrSampleData sF = feature->getData().createSampleData(sampleName);
                    if (feature->isInside(start, end, minOverlap)) {
                        sF->increaseReads();
                        if (end > feature->getEnd() && it2 != --(isoform->getFeatures().end()) && !countedIsoformBridge) {
                            sI->increaseBridgesReads();
                            countedIsoformBridge = true;
                        }
                        if (feature->getType().compare("exon") == 0 && !countedIsoformExon) {
                            sI->increaseExonReads();
                            countedIsoformExon = true;
                        }
                        if (feature->getType().compare("intron") == 0 && !countedIsoformIntron) {
                            sI->increaseIntronReads();
                            countedIsoformIntron = true;
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
    double TGene = 0.0;
    double TGeneExon = 0.0;
    double TGeneIntron = 0.0;
    double TGeneUnique = 0.0;
    double TGeneUniqueExon = 0.0;
    double TGeneUniqueIntron = 0.0;
    double TGeneFeature = 0.0;
    double TGeneUniqueFeature = 0.0;

    double TIsoform = 0.0;
    double TIsoformExon = 0.0;
    double TIsoformIntron = 0.0;
    double TIsoformFeature = 0.0;


    double TBridges = 0.0;


    for (auto cIt = genomeFactory.getChromosomes().begin(); cIt != genomeFactory.getChromosomes().end(); ++cIt) {
        c = cIt->second;
        for (auto gIt = c->getGenes().begin(); gIt != c->getGenes().end(); ++gIt) {
            g = *gIt;
            if (g->isProcessed()) {
                try {
                    SPtrSampleData s = g->getData().getSampleData(sampleName);
                    TGene += static_cast<double> (s->getReads()) / static_cast<double> (g->getLength());
                    if (g->getExonLength() > 0) {
                        TGeneExon += static_cast<double> (s->getExonReads()) / static_cast<double> (g->getExonLength());
                    }
                    if (g->getIntronLength() > 0) {
                        TGeneIntron += static_cast<double> (s->getIntronReads()) / static_cast<double> (g->getIntronLength());
                    }
                    if (g->getUniqueLength() > 0) {
                        TGeneUnique += static_cast<double> (s->getUniqueReads()) / static_cast<double> (g->getUniqueLength());
                    }
                    if (g->getUniqueExonLength() > 0) {
                        TGeneUniqueExon += static_cast<double> (s->getUniqueExonReads()) / static_cast<double> (g->getUniqueExonLength());
                    }
                    if (g->getUniqueIntronLength() > 0) {
                        TGeneUniqueIntron += static_cast<double> (s->getUniqueIntronReads()) / static_cast<double> (g->getUniqueIntronLength());
                    }

                    for (auto fIt = g->getFeatures().begin(); fIt != g->getFeatures().end(); ++fIt) {
                        f = *fIt;
                        SPtrSampleData sf = f->getData().getSampleData(sampleName);
                        TGeneFeature += static_cast<double> (sf->getReads()) / static_cast<double> (f->getLength());
                    }
                    for (auto fIt = g->getUniqueFeatures().begin(); fIt != g->getUniqueFeatures().end(); ++fIt) {
                        f = *fIt;
                        SPtrSampleData sf = f->getData().getSampleData(sampleName);
                        TGeneUniqueFeature += static_cast<double> (sf->getReads()) / static_cast<double> (f->getLength());
                    }
                } catch (exceptions::NotFoundException) {
                }

                for (auto iIt = g->getIsoforms().begin(); iIt != g->getIsoforms().end(); ++iIt) {
                    i = *iIt;
                    if (i->isProcessed()) {
                        try {
                            SPtrSampleData s = i->getData().createSampleData(sampleName);

                            TIsoform += static_cast<double> (s->getReads()) / static_cast<double> (i->getLength());
                            if (i->getExonLength() > 0) {
                                TIsoformExon += static_cast<double> (s->getExonReads()) / static_cast<double> (i->getExonLength());
                            }
                            if (i->getIntronLength() > 0) {
                                TIsoformIntron += static_cast<double> (s->getIntronReads()) / static_cast<double> (i->getIntronLength());
                            }
                            if (i->getFeatures().size() > 1) {
                                TBridges += static_cast<double> (s->getBridgesReads()) / static_cast<double> ((i->getFeatures().size() - 1));
                            }

                            for (auto fIt = i->getFeatures().begin(); fIt != i->getFeatures().end(); ++fIt) {
                                f = *fIt;
                                SPtrSampleData sf = f->getData().createSampleData(sampleName);
                                TIsoformFeature += static_cast<double> (sf->getReads()) / static_cast<double> (f->getLength());
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
                try {
                    SPtrSampleData s = g->getData().getSampleData(sampleName);
                    if (s->getReads() != 0) {
                        s->setTPM(static_cast<double> (s->getReads() * 1.0E6) / static_cast<double> (g->getLength() * TGene));
                    }
                    if (s->getExonReads() != 0) {
                        s->setExonTPM(static_cast<double> (s->getExonReads() * 1.0E6) / static_cast<double> (g->getExonLength() * TGeneExon));
                    }
                    if (s->getIntronReads() != 0) {
                        s->setIntronTPM(static_cast<double> (s->getIntronReads() * 1.0E6) / static_cast<double> (g->getIntronLength() * TGeneIntron));
                    }
                    if (s->getUniqueReads() != 0) {
                        s->setUniqueTPM(static_cast<double> (s->getUniqueReads() * 1.0E6) / static_cast<double> (g->getUniqueLength() * TGeneUnique));
                    }
                    if (s->getUniqueExonReads() != 0) {
                        s->setUniqueExonTPM(static_cast<double> (s->getUniqueExonReads() * 1.0E6) / static_cast<double> (g->getUniqueExonLength() * TGeneUniqueExon));
                    }
                    if (s->getUniqueIntronReads() != 0) {
                        s->setUniqueIntronTPM(static_cast<double> (s->getUniqueIntronReads() * 1.0E6) / static_cast<double> (g->getUniqueIntronLength() * TGeneUniqueIntron));
                    }

                    for (auto fIt = g->getFeatures().begin(); fIt != g->getFeatures().end(); ++fIt) {
                        f = *fIt;
                        try {
                            SPtrSampleData sf = f->getData().getSampleData(sampleName);
                            sf->setTPM(static_cast<double> (sf->getReads() * 1.0E6) / static_cast<double> (f->getLength() * TGeneFeature));
                        } catch (exceptions::NotFoundException) {
                        }
                    }
                    for (auto fIt = g->getUniqueFeatures().begin(); fIt != g->getUniqueFeatures().end(); ++fIt) {
                        f = *fIt;
                        try {
                            SPtrSampleData sf = f->getData().getSampleData(sampleName);
                            sf->setTPM(static_cast<double> (sf->getReads() * 1.0E6) / static_cast<double> (f->getLength() * TGeneUniqueFeature));
                        } catch (exceptions::NotFoundException) {
                        }
                    }
                } catch (exceptions::NotFoundException) {
                }

                for (auto iIt = g->getIsoforms().begin(); iIt != g->getIsoforms().end(); ++iIt) {
                    i = *iIt;
                    if (i->isProcessed()) {
                        try {
                            SPtrSampleData s = i->getData().getSampleData(sampleName);
                            if (s->getReads() != 0) {
                                s->setTPM(static_cast<double> (s->getReads() * 1.0E6) / static_cast<double> (i->getLength() * TIsoform));
                            }
                            if (s->getExonReads() != 0) {
                                s->setExonTPM(static_cast<double> (s->getExonReads() * 1.0E6) / static_cast<double> (i->getExonLength() * TIsoformExon));
                            }
                            if (s->getIntronReads() != 0) {
                                s->setIntronTPM(static_cast<double> (s->getIntronReads() * 1.0E6) / static_cast<double> (i->getIntronLength() * TIsoformIntron));
                            }

                            if (i->getFeatures().size() > 1) {
                                s->setBridgesTPM(static_cast<double> (s->getBridgesReads() * 1.0E6) / static_cast<double> ((i->getFeatures().size() - 1) * TBridges));
                            }
                            for (auto fIt = i->getFeatures().begin(); fIt != i->getFeatures().end(); ++fIt) {
                                f = *fIt;
                                try {
                                    SPtrSampleData sf = f->getData().createSampleData(sampleName);
                                    sf->setTPM(static_cast<double> (sf->getReads() * 1.0E6) / static_cast<double> (f->getLength() * TIsoformFeature));
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

int ReadFactory::processReadsFromBAM(std::string bamFileName, std::string sampleName, bool onlyProperlyPaired, uint16_t minMAPQ, uint16_t minOverlap) {
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
    fprintf(stderr, "\n");
    while (reader.GetNextAlignment(al)) {
        if (al.IsMapped() && al.MapQuality >= minMAPQ && !al.IsFailedQC()) {
            toRun = true;
            if (onlyProperlyPaired && !al.IsProperPair()) toRun = false;
            if (toRun) {
                std::set < std::pair<unsigned int, unsigned int>, coordinateLessCMP> read_coords;
                string chr = references[al.RefID].RefName;
                // al.Position position (0-based) where alignment starts (from Bamtools)
                start = end = al.Position;
                len = 0;
                start = al.Position;
                for (CigarOp c : al.CigarData) {
                    if (c.Type != 'S' && c.Type != 'H' && c.Type != 'P')
                        len += c.Length;
                    if (c.Type == 'N') {
                        if (end - start + 1 >= minOverlap) {
                            read_coords.insert(std::make_pair(start, end));
                        }
                        start = end + 1 + c.Length;
                        len = 0;
                    } else {
                        end = start + len - 1;
                    }
                }
                if (end - start + 1 >= minOverlap) {
                    read_coords.insert(std::make_pair(start, end));
                }
                if (!read_coords.empty()) {
                    processReadAtGenomeLevel(chr, sampleName, read_coords, minOverlap);
                    count++;
                }
            }
        }
    }
    reader.Close();

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

int ReadFactory::processBAMSAMFromDir(std::string dirName, bool onlyProperlyPaired, uint16_t minMAPQ, uint16_t minOverlap) {
    int count = 0;
    int totalCount = 0;
    struct dirent *dp;
    string BAMsufix(".bam");
    TimeUtils uTime;
    string sampleFileName;

    DIR *dirp = (DIR *) opendir(dirName.c_str());
    if (!dirp) {
        cerr << "Can't open directory: " << dirName << endl;
        exit(-1);
    }

    while ((dp = readdir(dirp)) != NULL) {
        string fName(dp->d_name);
        if (fName[0] != '.' && fName.size() >= 4) {
            if (fName.compare(fName.size() - BAMsufix.size(), BAMsufix.size(), BAMsufix) == 0) {
                string fileName = dirName + "/" + fName;
                sampleFileName = fName;
                sampleFileName.replace(fName.size() - BAMsufix.size(), BAMsufix.size(), "");
                uTime.setTime();
                cerr << "Processing sample: " << sampleFileName;
                samples.push_back(sampleFileName);
                cerr.flush();
                try {
                    count = processReadsFromBAM(fileName, sampleFileName, onlyProperlyPaired, minMAPQ, minOverlap);
                    totalCount += count;
                } catch (exceptions::NotFoundException) {
                    cerr << "Can't open file " << fileName << endl;
                    exit(-1);
                }
                cerr << " in " << uTime.getElapseTimeSec() << " seconds. Processed " << count << " reads" << endl;
            }
        }
    }
    closedir(dirp);

    return totalCount;
}

void ReadFactory::printResults(bool singleFile) {
    int count = 0;
    string sampleFileName;
    SPtrChromosomeNGS c;
    SPtrGeneNGS g;
    SPtrIsoformNGS i;
    SPtrFeatureNGS f;
    ofstream out_trans, ent_trans, out_gene, ent_gene, ent_gene_unique;

    for (auto sIt = samples.begin(); sIt != samples.end(); ++sIt) {
        sampleFileName = *sIt + "_transcripts.out";
        out_trans.open(sampleFileName);
        if (!out_trans.is_open()) {
            cerr << "Can't open file " << sampleFileName << endl;
            exit(-1);
        }
        out_trans << "Gene_Id\tTranscript_Id\tChr\tStart\tEnd\tLength\tReads\tTPM\tExonLength\tExonReads\tExonTPM\tIntronLength\tIntronReads\tIntronTPM" << endl;

        sampleFileName = *sIt + "_transcripts.ent";
        ent_trans.open(sampleFileName);
        ent_trans << "Gene_Id\tTranscript_Id\tChr\tType\tType_Number\tstart\tend\tLength\tReads\tTPM" << endl;
        if (!ent_trans.is_open()) {
            cerr << "Can't open file " << sampleFileName << endl;
            exit(-1);
        }

        sampleFileName = *sIt + "_genes.out";
        out_gene.open(sampleFileName);
        if (!out_gene.is_open()) {
            cerr << "Can't open file " << sampleFileName << endl;
            exit(-1);
        }
        out_gene << "Gene_Id\tChr\tStart\tEnd\tLength\tReads\tTPM\tExonLength\tExonReads\tExonTPM\tIntronLength\tIntronReads\tIntronTPM\tUniqueLegth\tUniqueReads\tUniqueTPM\tUniqueExonLength\tUniqueExonReads\tUniqueExonTPM\tUniqueIntronLength\tUniqueIntronReads\tUniqueIntronTPM" << endl;

        sampleFileName = *sIt + "_genes.ent";
        ent_gene.open(sampleFileName);
        ent_gene << "Gene_Id\tChr\tType\tType_Number\tstart\tend\tLength\tReads\tTPM" << endl;
        if (!ent_gene.is_open()) {
            cerr << "Can't open file " << sampleFileName << endl;
            exit(-1);
        }
        sampleFileName = *sIt + "_genes.uni";
        ent_gene_unique.open(sampleFileName);
        ent_gene_unique << "Gene_Id\tChr\tType\tType_Number\tstart\tend\tLength\tReads\tTPM" << endl;
        if (!ent_gene_unique.is_open()) {
            cerr << "Can't open file " << sampleFileName << endl;
            exit(-1);
        }

        for (auto cIt = genomeFactory.getChromosomes().begin(); cIt != genomeFactory.getChromosomes().end(); ++cIt) {
            c = cIt->second;
            for (auto it = c->getGenes().begin(); it != c->getGenes().end(); ++it) {
                g = *it;
                if (g->isProcessed()) {
                    out_gene << g->getId()
                            << "\t" << c->getId()
                            << "\t" << g->getStart()
                            << "\t" << g->getEnd()
                            << "\t" << g->getLength();
                    SPtrSampleData s = nullptr;
                    double reads = 0;
                    double tpm = 0.0;
                    double exonReads = 0;
                    double exonTPM = 0.0;
                    double intronReads = 0;
                    double intronTPM = 0.0;
                    double uniqueReads = 0;
                    double uniqueTPM = 0.0;
                    double uniqueExonReads = 0;
                    double uniqueExonTPM = 0.0;
                    double uniqueIntronReads = 0;
                    double uniqueIntronTPM = 0.0;

                    try {
                        s = g->getData().getSampleData(*sIt);
                        reads = s->getReads();
                        tpm = s->getTPM();
                        exonReads = s->getExonReads();
                        exonTPM = s->getExonTPM();
                        intronReads = s->getIntronReads();
                        intronTPM = s->getIntronTPM();
                        uniqueReads = s->getUniqueReads();
                        uniqueTPM = s->getUniqueTPM();
                        uniqueExonReads = s->getUniqueExonReads();
                        uniqueExonTPM = s->getUniqueExonTPM();
                        uniqueIntronReads = s->getUniqueIntronReads();
                        uniqueIntronTPM = s->getUniqueIntronTPM();
                    } catch (exceptions::NotFoundException) {
                    }

                    out_gene << "\t" << reads
                            << "\t" << tpm
                            << "\t" << g->getExonLength()
                            << "\t" << exonReads
                            << "\t" << exonTPM
                            << "\t" << g->getIntronLength()
                            << "\t" << intronReads
                            << "\t" << intronTPM
                            << "\t" << g->getUniqueLength()
                            << "\t" << uniqueReads
                            << "\t" << uniqueTPM
                            << "\t" << g->getUniqueExonLength()
                            << "\t" << uniqueExonReads
                            << "\t" << uniqueExonTPM
                            << "\t" << g->getUniqueIntronLength()
                            << "\t" << uniqueIntronReads
                            << "\t" << uniqueIntronTPM
                            << endl;

                    count = 1;
                    for (auto fIt = g->getFeatures().begin(); fIt != g->getFeatures().end(); ++fIt) {
                        f = *fIt;
                        try {
                            SPtrSampleData s = f->getData().getSampleData(*sIt);
                            ent_gene << g->getId()
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
                            ent_gene << g->getId()
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

                    count = 1;
                    for (auto fIt = g->getUniqueFeatures().begin(); fIt != g->getUniqueFeatures().end(); ++fIt) {
                        f = *fIt;
                        try {
                            SPtrSampleData s = f->getData().getSampleData(*sIt);
                            ent_gene_unique << g->getId()
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
                            ent_gene_unique << g->getId()
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

                    for (auto it1 = g->getIsoforms().begin(); it1 != g->getIsoforms().end(); ++it1) {
                        i = *it1;
                        if (i->isProcessed()) {
                            out_trans << g->getId()
                                    << "\t" << i->getId()
                                    << "\t" << c->getId()
                                    << "\t" << i->getStart()
                                    << "\t" << i->getEnd()
                                    << "\t" << i->getLength();
                            SPtrSampleData s;
                            double reads = 0;
                            double tpm = 0.0;
                            double exonReads = 0;
                            double exonTPM = 0.0;
                            double intronReads = 0;
                            double intronTPM = 0.0;
                            try {
                                SPtrSampleData s = i->getData().getSampleData(*sIt);
                                reads = s->getReads();
                                tpm = s->getTPM();
                                exonReads = s->getExonReads();
                                exonTPM = s->getExonTPM();
                                intronReads = s->getIntronReads();
                                intronTPM = s->getIntronTPM();
                            } catch (exceptions::NotFoundException) {
                            }
                            out_trans << "\t" << reads
                                    << "\t" << tpm
                                    << "\t" << i->getExonLength()
                                    << "\t" << exonReads
                                    << "\t" << exonTPM
                                    << "\t" << i->getIntronLength()
                                    << "\t" << intronReads
                                    << "\t" << intronTPM
                                    << endl;

                            count = 1;
                            for (auto fIt = i->getFeatures().begin(); fIt != i->getFeatures().end(); ++fIt) {
                                f = *fIt;
                                try {
                                    SPtrSampleData s = f->getData().getSampleData(*sIt);
                                    ent_trans << g->getId()
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
                                    ent_trans << g->getId()
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
        out_trans.close();
        ent_trans.close();
        out_gene.close();
        ent_gene.close();
        ent_gene_unique.close();
    }
}

void ReadFactory::createSIMSingleReadsIR(std::string outFileName, sequence::DNAContainer seqContainer, unsigned int numberFeat, unsigned int intronNumber, unsigned int len) {

    Random rng = Random();
    ofstream outputFile01(outFileName + "_genes.txt");
    ofstream outputFile02(outFileName + "_IR_genes.txt");
    ofstream outputFile1(outFileName + "_cond1.fa");
    ofstream outputFile2(outFileName + "_cond2.fa");
    for (ChromosomeUnMapItr<ReadData> cIt = genomeFactory.getChromosomes().begin();
            cIt != genomeFactory.getChromosomes().end(); ++cIt) {
        SPtrChromosomeNGS chromosome = cIt->second;

        if (chromosome->getId().size() <= 5) {
            try {
                SPtrDNA d = seqContainer.getDNAFromID(chromosome->getId());
                cout << "Getting " << chromosome->getId() << " chromosome sequence" << endl;
                cout << "Chromosome length " << d->getLength() << endl;
                set<pair<unsigned int, unsigned int>>genes_coord;
                std::vector < SPtrGene < ReadData>> genes;

                for (GeneMultiSetItr<ReadData> gIt = chromosome->getGenes().begin();
                        gIt != chromosome->getGenes().end(); ++gIt) {
                    SPtrGeneNGS gene = *gIt;
                    bool selected = false;
                    for (IsoformMultiSetItr<ReadData> iIt = gene->getIsoforms().begin(); iIt != gene->getIsoforms().end(); ++iIt) {
                        SPtrIsoformNGS isoform = *iIt;
                        if (isoform->getFeatures().size() >= (2 * numberFeat - 1)) {
                            selected = true;
                            break;
                        }
                    }
                    if (selected) {
                        selected = false;
                        for (auto genesItr = genes_coord.begin(); genesItr != genes_coord.end(); ++genesItr) {
                            pair<unsigned int, unsigned int> p = *genesItr;
                            if (gene->isInside(p.first, p.second, 0)) {
                                selected = true;
                                break;
                            }
                        }
                        if (!selected) {
                            if (gene->getFeatures().size() > 0) {
                                pair<unsigned int, unsigned int> coord = make_pair(gene->getStart(), gene->getEnd());
                                genes_coord.insert(coord);
                                genes.push_back(gene);
                                selected = false;
                                for (auto coordItr = gene->getFeatures().begin(); coordItr != gene->getFeatures().end(); ++coordItr) {
                                    SPtrFeatureNGS f = *coordItr;
                                    if (f->getLength() > len && f->getType() == "exon") {
                                        for (int k = 0; k < 2; k++) {
                                            int readNumber = f->getLength() * rng.DrawNumber(2, 4);
                                            for (int j = readNumber * 10000; j > 0; j--) {
                                                unsigned int min = f->getStart() - 7;
                                                unsigned int max = f->getEnd() + 8 - len;
                                                unsigned int start_pos = rng.DrawNumber(min, max);
                                                if (start_pos + len < d->getLength()) {
                                                    DNA read = d->newSegment(start_pos, len);
                                                    if (read.getSeq().find('N') == std::string::npos && read.getSeq().find('n') == std::string::npos) {
                                                        readNumber--;
                                                        selected = true;
                                                        if (k == 0) {
                                                            outputFile1 << ">" << read.getId() << ":" << start_pos << "-" << (start_pos + len - 1) << endl;
                                                        } else {
                                                            outputFile2 << ">" << read.getId() << ":" << start_pos << "-" << (start_pos + len - 1) << endl;
                                                        }
                                                        for (unsigned int i = 0; i < read.getLength(); i += 50) {
                                                            char t = 0;
                                                            if (i + 50 < read.getLength()) {
                                                                t = read.getSeq()[i + 50];
                                                                read.getSeq()[i + 50] = 0;
                                                            }
                                                            if (k == 0) {
                                                                outputFile1 << (read.getSeq().c_str() + i) << endl;
                                                            } else {
                                                                outputFile2 << (read.getSeq().c_str() + i) << endl;
                                                            }
                                                            if (i + 50 < read.getLength()) {
                                                                read.getSeq()[i + 50] = t;
                                                            }
                                                        }
                                                        if (readNumber == 0) break;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                if (selected) {
                                    outputFile01 << chromosome->getId() << "\t" << gene->getId() << endl;
                                }
                            }
                        }
                    }
                }

                if (!genes.empty()) {
                    /*
                     * Selecting genes for IR
                     */
                    cout << "There are " << genes.size() << " loaded" << endl;
                    cout << "20% of the genes (" << ((int) genes.size() / 5) << ") will get IR" << endl;
                    for (int ii = 0; ii < (int) genes.size() / 5;) {
                        unsigned int index = rng.DrawNumber(0, genes.size() - 1);
                        SPtrGeneNGS gene = genes[index];
                        if (gene->getFeatures().size() > intronNumber) {
                            unsigned int k = 0;
                            std::vector<unsigned int> index_selected;
                            for (int iterations = 0; iterations < 10000; iterations++) {
                                unsigned int index = rng.DrawNumber(0, gene->getFeatures().size() - 1);
                                if (std::find(index_selected.begin(), index_selected.end(), index) == index_selected.end()) {
                                    auto coordItr = gene->getFeatures().begin();
                                    std::advance(coordItr, index);
                                    SPtrFeatureNGS f = *coordItr;
                                    if (f->getLength() > len && f->getType() == "intron") {
                                        index_selected.push_back(index);
                                        int readNumber = f->getLength() * rng.DrawNumber(1, 2);
                                        for (int j = readNumber * 1000; j != 0; j--) {
                                            unsigned int min = f->getStart() - 7;
                                            unsigned int max = f->getEnd() + 8 - len;
                                            unsigned int start_pos = rng.DrawNumber(min, max);
                                            if (start_pos + len < d->getLength()) {
                                                DNA read = d->newSegment(start_pos, len);
                                                if (read.getSeq().find('N') == std::string::npos && read.getSeq().find('n') == std::string::npos) {
                                                    readNumber--;
                                                    outputFile2 << ">" << read.getId() << ":" << start_pos << "-" << (start_pos + len - 1) << endl;
                                                    for (unsigned int i = 0; i < read.getLength(); i += 50) {
                                                        char t = 0;
                                                        if (i + 50 < read.getLength()) {
                                                            t = read.getSeq()[i + 50];
                                                            read.getSeq()[i + 50] = 0;
                                                        }
                                                        outputFile2 << (read.getSeq().c_str() + i) << endl;
                                                        if (i + 50 < read.getLength()) {
                                                            read.getSeq()[i + 50] = t;
                                                        }
                                                    }
                                                    if (readNumber == 0) {
                                                        break;
                                                    }
                                                }
                                            }
                                        }
                                        k++;
                                        if (k == intronNumber || index_selected.size() == gene->getFeatures().size()) {
                                            outputFile02 << chromosome->getId() << "\t" << gene->getId() << endl;
                                            break;
                                        }
                                    }
                                }
                            }
                            ii++;
                        }
                    }
                }
            } catch (exceptions::NotFoundException ex) {
            }
        }

    }
    outputFile01.close();
    outputFile02.close();
    outputFile1.close();
    outputFile2.close();
}

void ReadFactory::loadTPMCalculatorGenesOutput(std::string dirName) {
    struct dirent *dp;
    string GENESentsufix("_genes.ent");
    string GENESoutsufix("_genes.out");
    string GENESunisufix("_genes.uni");
    TimeUtils uTime;
    string sampleFileName;
    SPtrGeneNGS g;

    DIR *dirp = (DIR *) opendir(dirName.c_str());
    if (!dirp) {
        cerr << "Can't open directory: " << dirName << endl;
        exit(-1);
    }

    while ((dp = readdir(dirp)) != NULL) {
        string fName(dp->d_name);
        if (fName[0] != '.' && fName.size() >= 10) {
            if (fName.compare(fName.size() - GENESentsufix.size(), GENESentsufix.size(), GENESentsufix) == 0) {
                string fileName = dirName + "/" + fName;
                sampleFileName = fName;
                sampleFileName.replace(fName.size() - GENESentsufix.size(), GENESentsufix.size(), "");
                uTime.setTime();
                cerr << "Processing sample ent file: " << sampleFileName;
                samples.push_back(sampleFileName);
                cerr.flush();
                try {
                    parsers::TextParser fParser;
                    fParser.setFileToParse(fileName);
                    while (fParser.iterate("#", "\t")) {
                        if (fParser.getWords().size() != 9) {
                            std::cout << "\nGenes TPM ent output file with wrong number of fields. It should be 9 tab separated fields" << std::endl;
                            std::cout << "Words size: " << fParser.getWords().size() << std::endl;
                            std::cout << "Line size: " << fParser.getLine().size() << std::endl;
                            std::cout << fParser.getLine() << std::endl;
                            exit(-1);
                        }
                        if (fParser.getWords()[0] != "Gene_Id") {
                            try {
                                g = genomeFactory.findGene(fParser.getWords()[1], fParser.getWords()[0]);
                                SPtrSampleData s = g->getData().createSampleData(sampleFileName);
                                g->setProcessed(true);
                                int uNumber = atoi(fParser.getWords()[3].c_str()) - 1;
                                int uIndex = 0;
                                SPtrFeatureNGS f;
                                for (auto fIt : g->getFeatures()) {
                                    if (uIndex == uNumber) {
                                        f = fIt;
                                        break;
                                    }
                                    uIndex++;
                                }
                                if (f) {
                                    SPtrSampleData sf = f->getData().createSampleData(sampleFileName);
                                    sf->increaseReads(atoi(fParser.getWords()[7].c_str()));
                                    sf->setTPM(atof(fParser.getWords()[8].c_str()));
                                } else {
                                    std::cerr << "Feature not inserted" << std::endl;
                                    exit(-1);
                                }
                            } catch (exceptions::NotFoundException ex) {
                                std::cerr << "Chromosome: " << fParser.getWords()[1] << " not include in the GTF" << std::endl;
                                exit(-1);
                            }
                        }
                    }
                } catch (exceptions::FileHandledException ex) {
                    std::cerr << "Error parsing file: " << fileName << std::endl;
                    exit(-1);
                }
                cerr << " in " << uTime.getElapseTimeSec() << " seconds." << endl;
            }
            if (fName.compare(fName.size() - GENESunisufix.size(), GENESunisufix.size(), GENESunisufix) == 0) {
                string fileName = dirName + "/" + fName;
                sampleFileName = fName;
                sampleFileName.replace(fName.size() - GENESentsufix.size(), GENESentsufix.size(), "");
                uTime.setTime();
                cerr << "Processing sample uni file: " << sampleFileName;
                samples.push_back(sampleFileName);
                cerr.flush();
                try {
                    parsers::TextParser fParser;
                    fParser.setFileToParse(fileName);
                    while (fParser.iterate("#", "\t")) {
                        if (fParser.getWords().size() != 9) {
                            std::cout << "\nGenes TPM ent output file with wrong number of fields. It should be 9 tab separated fields" << std::endl;
                            std::cout << "Words size: " << fParser.getWords().size() << std::endl;
                            std::cout << "Line size: " << fParser.getLine().size() << std::endl;
                            std::cout << fParser.getLine() << std::endl;
                            exit(-1);
                        }
                        if (fParser.getWords()[0] != "Gene_Id") {
                            try {
                                g = genomeFactory.findGene(fParser.getWords()[1], fParser.getWords()[0]);
                                SPtrSampleData s = g->getData().createSampleData(sampleFileName);
                                g->setProcessed(true);
                                int uNumber = atoi(fParser.getWords()[3].c_str()) - 1;
                                int uIndex = 0;
                                SPtrFeatureNGS f;
                                for (auto fIt : g->getUniqueFeatures()) {
                                    if (uIndex == uNumber) {
                                        f = fIt;
                                        break;
                                    }
                                    uIndex++;
                                }
                                if (f) {
                                    SPtrSampleData sf = f->getData().createSampleData(sampleFileName);
                                    sf->increaseReads(atoi(fParser.getWords()[7].c_str()));
                                    sf->setTPM(atof(fParser.getWords()[8].c_str()));
                                } else {
                                    std::cerr << "Feature not inserted" << std::endl;
                                    exit(-1);
                                }
                            } catch (exceptions::NotFoundException ex) {
                                std::cerr << "Chromosome: " << fParser.getWords()[1] << " not include in the GTF" << std::endl;
                                exit(-1);
                            }
                        }
                    }
                } catch (exceptions::FileHandledException ex) {
                    std::cerr << "Error parsing file: " << fileName << std::endl;
                    exit(-1);
                }
                cerr << " in " << uTime.getElapseTimeSec() << " seconds." << endl;
            }
        }
        if (fName[0] != '.' && fName.size() >= 10) {
            if (fName.compare(fName.size() - GENESoutsufix.size(), GENESoutsufix.size(), GENESoutsufix) == 0) {
                string fileName = dirName + "/" + fName;
                sampleFileName = fName;
                sampleFileName.replace(fName.size() - GENESoutsufix.size(), GENESoutsufix.size(), "");
                uTime.setTime();
                cerr << "Processing sample out file: " << sampleFileName;
                cerr.flush();
                try {
                    parsers::TextParser fParser;
                    fParser.setFileToParse(fileName);
                    while (fParser.iterate("#", "\t")) {
                        if (fParser.getWords().size() != 20) {
                            std::cerr << "Genes TPM out output file with wrong number of fields. It should be 17 tab separated fields or 11 for the old version" << std::endl;
                            exit(-1);
                        }
                        if (fParser.getWords()[0] != "Gene_Id") {
                            try {
                                g = genomeFactory.findGene(fParser.getWords()[1], fParser.getWords()[0]);
                                SPtrSampleData s = g->getData().createSampleData(sampleFileName);
                                g->setProcessed(true);
                                s->increaseReads(atoi(fParser.getWords()[5].c_str()));
                                s->setTPM(atof(fParser.getWords()[6].c_str()));

                                s->increaseExonReads(atoi(fParser.getWords()[8].c_str()));
                                s->setExonTPM(atof(fParser.getWords()[9].c_str()));

                                s->increaseIntronReads(atoi(fParser.getWords()[11].c_str()));
                                s->setIntronTPM(atof(fParser.getWords()[12].c_str()));

                                s->increaseUniqueReads(atoi(fParser.getWords()[14].c_str()));
                                s->setUniqueTPM(atof(fParser.getWords()[15].c_str()));

                                s->increaseUniqueExonReads(atoi(fParser.getWords()[17].c_str()));
                                s->setUniqueExonTPM(atof(fParser.getWords()[18].c_str()));

                                s->increaseUniqueIntronReads(atoi(fParser.getWords()[20].c_str()));
                                s->setUniqueIntronTPM(atof(fParser.getWords()[21].c_str()));
                            } catch (exceptions::NotFoundException ex) {
                                std::cerr << "Chromosome: " << fParser.getWords()[1] << " not include in the GTF" << std::endl;
                                exit(-1);
                            }
                        }
                    }
                } catch (exceptions::FileHandledException ex) {
                    std::cerr << "Error parsing file: " << fileName << std::endl;
                    exit(-1);
                }
                cerr << " in " << uTime.getElapseTimeSec() << " seconds." << endl;
            }
        }
    }
    closedir(dirp);
}

void ReadFactory::printResultsMatrix(std::string output_name, std::vector<std::string> tpmColumns) {
    SPtrChromosomeNGS c;
    SPtrGeneNGS g;
    SPtrIsoformNGS i;
    SPtrFeatureNGS f;

    for (auto col : tpmColumns) {
        ofstream out;
        out.open(output_name + "_" + col + ".txt");
        if (!out.is_open()) {
            cerr << "Can't open file " << output_name << endl;
            exit(-1);
        }
        out << "Gene_Id";
        for (auto sIt = samples.begin(); sIt != samples.end(); ++sIt) {
            out << "\t" << *sIt;
        }
        out << endl;

        for (auto cIt = genomeFactory.getChromosomes().begin(); cIt != genomeFactory.getChromosomes().end(); ++cIt) {
            c = cIt->second;
            for (auto it = c->getGenes().begin(); it != c->getGenes().end(); ++it) {
                g = *it;
                if (g->isProcessed()) {
                    out << g->getId();
                    for (auto sIt = samples.begin(); sIt != samples.end(); ++sIt) {
                        try {
                            SPtrSampleData s = g->getData().getSampleData(*sIt);
                            out << "\t" << s->getValueFromColumn(col);
                        } catch (exceptions::NotFoundException) {
                            out << "\t0";
                        }
                    }
                    out << endl;
                }
            }
        }
        out.close();
    }


}
