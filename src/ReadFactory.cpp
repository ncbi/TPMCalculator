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

void ReadFactory::processReadAtGenomeLevel(std::string chrName, std::string sampleName, unsigned int start, unsigned int end) {
    bool done = false;
    GeneMultiSetNGS::iterator geneIt;
    try {
        geneIt = genomeFactory.findGeneUpperBound(chrName, start, end);
        for (auto it = geneIt;; --it) {
            this->processReadAtGeneLevel(*it, sampleName, start, end);
            if (done) break;
            if (std::distance(it, geneIt) > 6 && (*it)->getEnd() < start) done = true;
            if (it == genomeFactory.getCurrentChr()->getGenes().begin()) break;
        }
    } catch (exceptions::NotFoundException) {
    }
}

void ReadFactory::processReadAtGenomeLevelUnique(std::string chrName, std::string sampleName, unsigned int start, unsigned int end) {
    bool done = false;
    GeneMultiSetNGS::iterator geneIt;
    try {
        geneIt = genomeFactory.findGeneUpperBound(chrName, start, end);
        for (auto it = geneIt;; --it) {
            this->processReadAtGeneLevelUnique(*it, sampleName, start, end);
            if (done) break;
            if (std::distance(it, geneIt) > 6 && (*it)->getEnd() < start) done = true;
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

void ReadFactory::processReadAtGeneLevelUnique(SPtrGeneNGS gene, std::string sampleName, unsigned int start, unsigned int end) {
    if (!gene->isInside(start, end, 8)) return;
    if (!gene->isProcessed()) gene->setProcessed(true);
    SPtrSampleData s = gene->getData().createSampleData(sampleName);
    for (auto it = gene->getUniquefeatures().begin(); it != gene->getUniquefeatures().end(); ++it) {
        if ((*it)->isInside(start, end, 8)) {
            (*it)->getData().increaseReads(sampleName);
            if (end > (*it)->getEnd() && it != --(gene->getUniquefeatures().end())) {
                s->increaseBridgeReads();
            }
        }
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
    SPtrFeatureNGS f;

    for (auto cIt = genomeFactory.getChromosomes().begin(); cIt != genomeFactory.getChromosomes().end(); ++cIt) {
        c = cIt->second;
        for (auto gIt = c->getGenes().begin(); gIt != c->getGenes().end(); ++gIt) {
            g = *gIt;
            if (g->isProcessed()) {
                try {
                    SPtrSampleData s = g->getData().createSampleData(sampleName);
                    for (auto fIt = g->getUniquefeatures().begin(); fIt != g->getUniquefeatures().end(); ++fIt) {
                        f = *fIt;
                        try {
                            if (f->getData().getSampleData(sampleName)->getReads() != 0 && f->getType().compare("exon") == 0) {
                                s->increaseReads(f->getData().getSampleData(sampleName)->getReads());
                            }
                        } catch (exceptions::NotFoundException) {
                        }
                    }
                } catch (exceptions::NotFoundException) {
                }
            } else {
                g->setProcessed(true);
                try {
                    SPtrSampleData s = g->getData().createSampleData(sampleName);
                    for (auto fIt = g->getUniquefeatures().begin(); fIt != g->getUniquefeatures().end(); ++fIt) {
                        f = *fIt;
                        try {
                            SPtrSampleData sf = f->getData().createSampleData(sampleName);
                        } catch (exceptions::NotFoundException) {
                        }
                    }
                } catch (exceptions::NotFoundException) {
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
    double TIsoform = 0.0;
    double TGeneFeature = 0.0;
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

                    for (auto fIt = g->getUniquefeatures().begin(); fIt != g->getUniquefeatures().end(); ++fIt) {
                        f = *fIt;
                        try {
                            SPtrSampleData sf = f->getData().getSampleData(sampleName);
                            TGeneFeature += static_cast<double> (sf->getReads()) / static_cast<double> (f->getLength());
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
                } catch (exceptions::NotFoundException) {
                }

                for (auto iIt = g->getIsoforms().begin(); iIt != g->getIsoforms().end(); ++iIt) {
                    i = *iIt;
                    if (i->isProcessed()) {

                        try {
                            SPtrSampleData s = i->getData().createSampleData(sampleName);

                            TIsoform += static_cast<double> (s->getReads()) / static_cast<double> (i->getLength());
                            if (i->getFeatures().size() > 1) {
                                TBridges += static_cast<double> (s->getBridgeReads()) / static_cast<double> ((i->getFeatures().size() - 1));
                            }
                            for (auto fIt = i->getFeatures().begin(); fIt != i->getFeatures().end(); ++fIt) {
                                f = *fIt;
                                try {
                                    SPtrSampleData sf = f->getData().createSampleData(sampleName);
                                    TIsoformFeature += static_cast<double> (sf->getReads()) / static_cast<double> (f->getLength());
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
                    s->setTPM(static_cast<double> (s->getReads() * 1.0E6) / static_cast<double> (i->getLength() * TGene));

                    if (s->getExonLength() != 0) {
                        s->setTPMExon(static_cast<double> (s->getExonReads() * 1.0E6) / static_cast<double> (s->getExonLength() * TGeneFeature));
                    }
                    if (s->getIntronLength() != 0) {
                        s->setTPMIntron(static_cast<double> (s->getIntronReads() * 1.0E6) / static_cast<double> (s->getIntronLength() * TGeneFeature));
                    }

                    for (auto fIt = g->getUniquefeatures().begin(); fIt != g->getUniquefeatures().end(); ++fIt) {
                        f = *fIt;
                        try {
                            SPtrSampleData sf = f->getData().getSampleData(sampleName);
                            sf->setTPM(static_cast<double> (sf->getReads() * 1.0E6) / static_cast<double> (f->getLength() * TGeneFeature));
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
                            s->setTPM(static_cast<double> (s->getReads() * 1.0E6) / static_cast<double> (i->getLength() * TIsoform));

                            if (i->getFeatures().size() > 1) {
                                s->setTPMBridges(static_cast<double> (s->getBridgeReads() * 1.0E6) / static_cast<double> ((i->getFeatures().size() - 1) * TBridges));
                            }

                            if (s->getExonLength() != 0) {
                                s->setTPMExon(static_cast<double> (s->getExonReads() * 1.0E6) / static_cast<double> (s->getExonLength() * TIsoformFeature));
                            }
                            if (s->getIntronLength() != 0) {
                                s->setTPMIntron(static_cast<double> (s->getIntronReads() * 1.0E6) / static_cast<double> (s->getIntronLength() * TIsoformFeature));
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

int ReadFactory::processReadsFromBAM(std::string bamFileName, std::string sampleName, bool onlyProperlyPaired) {
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
        toRun = false;
        if (al.IsMapped()) toRun = true;
        if (toRun && onlyProperlyPaired && !al.IsProperPair()) toRun = false;
        if (toRun) {
            string chr = references[al.RefID].RefName;
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
                            processReadAtGenomeLevel(chr, sampleName, start, end);
                            processReadAtGenomeLevelUnique(chr, sampleName, start, end);
                        }
                        start = end + 1 + c.Length;
                        len = 0;
                    } else {
                        end = start + len - 1;
                    }
                }
                if (len >= 8) {
                    processReadAtGenomeLevel(chr, sampleName, start, end);
                    processReadAtGenomeLevelUnique(chr, sampleName, start, end);
                }
            } else {
                processReadAtGenomeLevel(chr, sampleName, al.Position, al.GetEndPosition(true, true));
                processReadAtGenomeLevelUnique(chr, sampleName, al.Position, al.GetEndPosition(true, true));
            }
            count++;
        }
    }

    reader.Close();

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

int ReadFactory::processBAMSAMFromDir(std::string dirName, bool onlyProperlyPaired) {
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
                    count = processReadsFromBAM(fileName, sampleFileName, onlyProperlyPaired);
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
    ofstream out, ent, out_unique, ent_unique, out_gene;

    for (auto sIt = samples.begin(); sIt != samples.end(); ++sIt) {
        sampleFileName = *sIt + "_transcripts.out";
        out.open(sampleFileName);
        if (!out.is_open()) {
            cerr << "Can't open file " << sampleFileName << endl;
            exit(-1);
        }
        out << "Gene_Id\tTranscript_Id\tChr\tLength\tCount_Reads\tTPM\tExon_Length\tExon_Count_Reads\tExon_TPM\tIntron_Length\tIntron_Count_Reads\tIntron_TPM" << endl;

        sampleFileName = *sIt + "_transcripts.ent";
        ent.open(sampleFileName);
        ent << "Gene_Id\tTranscript_Id\tChr\tType\tType_Number\tstart\tend\tLength\tCount_Reads\tTPM" << endl;
        if (!ent.is_open()) {
            cerr << "Can't open file " << sampleFileName << endl;
            exit(-1);
        }

        sampleFileName = *sIt + "_genes.out";
        out_unique.open(sampleFileName);
        if (!out_unique.is_open()) {
            cerr << "Can't open file " << sampleFileName << endl;
            exit(-1);
        }
        out_unique << "Gene_Id\tChr\tLength\tCount_Reads\tTPM\tExon_Length\tExon_Count_Reads\tExon_TPM\tIntron_Length\tIntron_Count_Reads\tIntron_TPM" << endl;

        sampleFileName = *sIt + "_genes.ent";
        ent_unique.open(sampleFileName);
        ent_unique << "Gene_Id\tChr\tType\tType_Number\tstart\tend\tLength\tCount_Reads\tTPM" << endl;
        if (!ent_unique.is_open()) {
            cerr << "Can't open file " << sampleFileName << endl;
            exit(-1);
        }

        for (auto cIt = genomeFactory.getChromosomes().begin(); cIt != genomeFactory.getChromosomes().end(); ++cIt) {
            c = cIt->second;
            for (auto it = c->getGenes().begin(); it != c->getGenes().end(); ++it) {
                g = *it;
                if (g->isProcessed()) {
                    out_unique << g->getId()
                            << "\t" << c->getId()
                            << "\t" << g->getLength();
                    try {
                        SPtrSampleData s = g->getData().getSampleData(*sIt);
                        out_unique << "\t" << s->getReads()
                                << "\t" << s->getTPM()
                                << "\t" << s->getExonLength()
                                << "\t" << s->getExonReads()
                                << "\t" << s->getTPMExon()
                                << "\t" << s->getIntronLength()
                                << "\t" << s->getIntronReads()
                                << "\t" << s->getTPMIntron()
                                << endl;
                    } catch (exceptions::NotFoundException) {
                        out_unique << "\t0\t0.000000\t0\t0\t0.000000\t0\t0\t0.000000" << endl;
                    }
                    count = 1;
                    for (auto fIt = g->getUniquefeatures().begin(); fIt != g->getUniquefeatures().end(); ++fIt) {
                        f = *fIt;
                        try {
                            SPtrSampleData s = f->getData().getSampleData(*sIt);
                            ent_unique << g->getId()
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
                            ent_unique << g->getId()
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
        out_unique.close();
        ent_unique.close();
    }

    if (!singleFile) {
        sampleFileName = "transcripts_data_per_samples.txt";
        out.open(sampleFileName);
        if (!out.is_open()) {
            cerr << "Can't open file transcripts_data_per_samples.txt" << endl;
            exit(-1);
        }

        sampleFileName = "genes_data_per_samples.txt";
        out_gene.open(sampleFileName);
        if (!out_gene.is_open()) {
            cerr << "Can't open file genes_data_per_samples.txt" << endl;
            exit(-1);
        }

        out << "Gene_Id\tTranscript_Id";
        out_gene << "Gene_Id";
        for (auto sIt = samples.begin(); sIt != samples.end(); ++sIt) {
            out << "\t" << *sIt << "_exon_count\t" << *sIt << "_exon_TPM";
            out << "\t" << *sIt << "_intron_count\t" << *sIt << "_intron_TPM";
            out_gene << "\t" << *sIt << "_exon_count\t" << *sIt << "_exon_TPM";
            out_gene << "\t" << *sIt << "_intron_count\t" << *sIt << "_intron_TPM";
        }
        out << endl;
        out_gene << endl;

        for (auto cIt = genomeFactory.getChromosomes().begin(); cIt != genomeFactory.getChromosomes().end(); ++cIt) {
            c = cIt->second;
            for (auto it = c->getGenes().begin(); it != c->getGenes().end(); ++it) {
                g = *it;
                if (g->isProcessed()) {
                    out_gene << g->getId();
                    for (auto sIt = samples.begin(); sIt != samples.end(); ++sIt) {
                        try {
                            SPtrSampleData s = g->getData().getSampleData(*sIt);
                            out_gene << "\t" << s->getExonReads() << "\t" << s->getTPMExon();
                            out_gene << "\t" << s->getIntronReads() << "\t" << s->getTPMIntron();
                        } catch (exceptions::NotFoundException) {
                            out_gene << "\t0\t0.000000";
                            out_gene << "\t0\t0.000000";
                        }
                    }
                    out_gene << endl;

                    for (auto it1 = g->getIsoforms().begin(); it1 != g->getIsoforms().end(); ++it1) {
                        i = *it1;
                        if (i->isProcessed()) {
                            out << g->getId() << "\t" << i->getId();
                            for (auto sIt = samples.begin(); sIt != samples.end(); ++sIt) {
                                try {
                                    SPtrSampleData s = i->getData().getSampleData(*sIt);
                                    out << "\t" << s->getExonReads() << "\t" << s->getTPMExon();
                                    out << "\t" << s->getIntronReads() << "\t" << s->getTPMIntron();
                                } catch (exceptions::NotFoundException) {
                                    out << "\t0\t0.000000";
                                    out << "\t0\t0.000000";
                                }
                            }
                            out << endl;
                        }
                    }
                }
            }
        }
        out.close();
        out_gene.close();
    }
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
                            if (gene->getUniquefeatures().size() > 0) {
                                pair<unsigned int, unsigned int> coord = make_pair(gene->getStart(), gene->getEnd());
                                genes_coord.insert(coord);
                                genes.push_back(gene);
                                selected = false;
                                for (auto coordItr = gene->getUniquefeatures().begin(); coordItr != gene->getUniquefeatures().end(); ++coordItr) {
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
                        if (gene->getUniquefeatures().size() > intronNumber) {
                            unsigned int k = 0;
                            std::vector<unsigned int> index_selected;
                            for (int iterations = 0; iterations < 10000; iterations++) {
                                unsigned int index = rng.DrawNumber(0, gene->getUniquefeatures().size() - 1);
                                if (std::find(index_selected.begin(), index_selected.end(), index) == index_selected.end()) {
                                    auto coordItr = gene->getUniquefeatures().begin();
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
                                        if (k == intronNumber || index_selected.size() == gene->getUniquefeatures().size()) {
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
                cout << ex.what() << endl;
            }
        }

    }
    outputFile01.close();
    outputFile02.close();
    outputFile1.close();
    outputFile2.close();
}

void ReadFactory::loadTPMCalculatorGenesOutput(std::string dirName) {
    int count = 0;
    struct dirent *dp;
    string GENESentsufix("_genes.ent");
    string GENESoutsufix("_genes.out");
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
                                for (auto fIt : g->getUniquefeatures()) {
                                    if (uIndex == uNumber) {
                                        f = fIt;
                                        break;
                                    }
                                    uIndex++;
                                }
                                if (f) {
                                    //                                SPtrFeature<ReadData> f = *std::next(g->getUniquefeatures().begin(), atoi(fParser.getWords()[3].c_str()) - 1);
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
                        if (fParser.getWords().size() != 11) {
                            std::cerr << "Genes TPM out output file with wrong number of fields. It should be 11 tab separated fields" << std::endl;
                            exit(-1);
                        }
                        if (fParser.getWords()[0] != "Gene_Id") {
                            try {
                                g = genomeFactory.findGene(fParser.getWords()[1], fParser.getWords()[0]);
                                SPtrSampleData s = g->getData().createSampleData(sampleFileName);
                                g->setProcessed(true);
                                s->increaseReads(atoi(fParser.getWords()[3].c_str()));
                                s->setTPM(atof(fParser.getWords()[4].c_str()));
                                s->increaseExonLength(atoi(fParser.getWords()[5].c_str()));
                                s->increaseExonReads(atoi(fParser.getWords()[6].c_str()));
                                s->setTPMExon(atof(fParser.getWords()[7].c_str()));
                                s->increaseIntronLength(atoi(fParser.getWords()[8].c_str()));
                                s->increaseIntronReads(atoi(fParser.getWords()[9].c_str()));
                                s->setTPMIntron(atof(fParser.getWords()[10].c_str()));
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
