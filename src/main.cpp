#include <iostream>
#include <ctime>
#include <cstdint>
#include <set>
#include <unordered_map>
#include <algorithm>
#include <memory>
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

using namespace std;
using namespace ngs;
using namespace sequence;
using namespace genome;
using namespace BamTools;

Global *Global::s_instance = 0;

void print_usage(char *program_name, int exit_code) {
    cerr << "\n********************************************************************************\n";
    cerr << "\nUsage: " << program_name;
    cerr << "\n\n" << program_name << " options:\n\n";
    cerr << "-v    Print info\n";
    cerr << "-h    Display this usage information.\n";
    cerr << "-g    GTF file\n";
    cerr << "-d    Directory with the BAM files\n";
    cerr << "-b    BAM file\n";
    cerr << "-k    Gene key to use from GTF file. Default: gene_id\n";
    cerr << "-t    Transcript key to use from GTF file. Default: transcript_id\n";
    cerr << "-c    Smaller size allowed for an intron created for genes. Default: 16. We recommend to use the reads length\n";
    cerr << "-p    Use only properly paired reads. Default: No. Recommended for paired-end reads.\n";
    cerr << "-q    Minimum MAPQ value to filter out reads. Default: 0. This value depends on the aligner MAPQ value.\n";
    cerr << "-o    Minimum overlap between a reads and a feature. Default: 8.\n";
    cerr << "-e    Extended output. This will include transcript level TPM values. Default: No.\n";
    cerr << "-a    Print out all features with read counts equal to zero. Default: No.\n";
    cerr << "\n********************************************************************************\n";
    cerr << "\n                        Roberto Vera Alvarez, PhD\n";
    cerr << "                      Emails: veraalva@ncbi.nlm.nih.gov\n\n";
    cerr << "********************************************************************************\n";
    exit(exit_code);
}

int main(int argc, char *argv[]) {
    int count = 0;
    TimeUtils uTime;
    string gtfFileName;
    string bamDirName;
    string bamFileName;
    string geneNameKey = "gene_id";
    string transcriptNameKey = "transcript_id";
    int intronCutOff = 16;
    uint16_t minMAPQ = 0;
    uint16_t minOverlap = 8;
    bool onlyProperlyPaired = false;
    bool singleFile = false;
    bool extendedOutput = false;
    bool all_feat = false;
    set<string>features = {"exon"};
    unordered_map<string, string> featuresToCreate = {
        {"exon", "intron"}
    };
    ReadFactory readFactory;

    if (argc == 1) {
        print_usage(argv[0], 0);
    }

    for (int i = 1; i < argc; i++) {
        string option(argv[i]);
        if (option.compare(0, 1, "-") == 0 && option.compare(1, 1, "-") != 0 && option.size() == 2) {
            if (option.compare(1, 1, "h") == 0) {
                print_usage(argv[0], 0);
            } else if (option.compare(1, 1, "v") == 0) {
                Global::instance()->setVerbose(1);
            } else if (option.compare(1, 1, "g") == 0) {
                i++;
                if (i < argc) {
                    gtfFileName = argv[i];
                    if (gtfFileName.compare(0, 1, "-") == 0) {
                        cerr << "Option g require an argument" << endl;
                        print_usage(argv[0], -1);
                    }
                } else {
                    cerr << "Option g require an argument" << endl;
                    print_usage(argv[0], -1);
                }
            } else if (option.compare(1, 1, "d") == 0) {
                i++;
                if (i < argc) {
                    bamDirName = argv[i];
                    if (bamDirName.compare(0, 1, "-") == 0) {
                        cerr << "Option d require an argument" << endl;
                        print_usage(argv[0], -1);
                    }
                } else {
                    cerr << "Option d require an argument" << endl;
                    print_usage(argv[0], -1);
                }
            } else if (option.compare(1, 1, "b") == 0) {
                i++;
                if (i < argc) {
                    bamFileName = argv[i];
                    if (bamFileName.compare(0, 1, "-") == 0) {
                        cerr << "Option b require an argument" << endl;
                        print_usage(argv[0], -1);
                    }
                } else {
                    cerr << "Option b require an argument" << endl;
                    print_usage(argv[0], -1);
                }
            } else if (option.compare(1, 1, "k") == 0) {
                i++;
                if (i < argc) {
                    geneNameKey = argv[i];
                    if (geneNameKey.compare(0, 1, "-") == 0) {
                        cerr << "Option k require an argument" << endl;
                        print_usage(argv[0], -1);
                    }
                } else {
                    cerr << "Option k require an argument" << endl;
                    print_usage(argv[0], -1);
                }
            } else if (option.compare(1, 1, "t") == 0) {
                i++;
                if (i < argc) {
                    transcriptNameKey = argv[i];
                    if (transcriptNameKey.compare(0, 1, "-") == 0) {
                        cerr << "Option t require an argument" << endl;
                        print_usage(argv[0], -1);
                    }
                } else {
                    cerr << "Option t require an argument" << endl;
                    print_usage(argv[0], -1);
                }
            } else if (option.compare(1, 1, "c") == 0) {
                i++;
                if (i < argc) {
                    string argument(argv[i]);
                    if (argument.compare(0, 1, "-") == 0) {
                        cerr << "Option c require an argument" << endl;
                        print_usage(argv[0], -1);
                    }
                    intronCutOff = atoi(argv[i]);
                } else {
                    cerr << "Option c require an argument" << endl;
                    print_usage(argv[0], -1);
                }
            } else if (option.compare(1, 1, "q") == 0) {
                i++;
                if (i < argc) {
                    string argument(argv[i]);
                    if (argument.compare(0, 1, "-") == 0) {
                        cerr << "Option q require an argument" << endl;
                        print_usage(argv[0], -1);
                    }
                    minMAPQ = static_cast<uint16_t> (atoi(argv[i]));
                } else {
                    cerr << "Option q require an argument" << endl;
                    print_usage(argv[0], -1);
                }
            } else if (option.compare(1, 1, "o") == 0) {
                i++;
                if (i < argc) {
                    string argument(argv[i]);
                    if (argument.compare(0, 1, "-") == 0) {
                        cerr << "Option o require an argument" << endl;
                        print_usage(argv[0], -1);
                    }
                    minOverlap = static_cast<uint16_t> (atoi(argv[i]));
                } else {
                    cerr << "Option o require an argument" << endl;
                    print_usage(argv[0], -1);
                }
            } else if (option.compare(1, 1, "p") == 0) {
                onlyProperlyPaired = true;
            } else if (option.compare(1, 1, "e") == 0) {
                extendedOutput = true;
            } else if (option.compare(1, 1, "a") == 0) {
                all_feat = true;
            } else {
                cerr << "Unsupported option: " << option << endl;
                print_usage(argv[0], -1);
            }
        } else {
            cerr << "Unsupported option: " << option << endl;
            print_usage(argv[0], -1);
        }
    }

    if (gtfFileName.empty()) {
        cerr << "\nGTF is required. See -g option" << endl;
        print_usage(argv[0], -1);
    }

    if (bamDirName.empty() && bamFileName.empty()) {
        cerr << "\nDirectory with the BAM files or a BAM file is required. See -d or -b options" << endl;
        print_usage(argv[0], -1);
    }

    uTime.setTime();
    cerr << "Reading GTF file ... " << endl;
    readFactory.getGenomeFactory().setIntronCutOff(intronCutOff);
    readFactory.getGenomeFactory().processGTFFile(gtfFileName, geneNameKey, transcriptNameKey, features, featuresToCreate);
    cerr << "Done in " << uTime.getElapseTimeSec() << " seconds" << endl;

    if (!bamDirName.empty()) {
        uTime.setTime();
        cerr << "Parsing BAM files" << endl;
        fflush(NULL);
        count = readFactory.processBAMSAMFromDir(bamDirName, onlyProperlyPaired, minMAPQ, minOverlap);
        cerr << count << " reads processed in " << uTime.getElapseTimeSec() << " seconds" << endl;
        fflush(NULL);
    } else if (!bamFileName.empty()) {
        singleFile = true;
        string fileName = bamFileName;
        string sampleName;
        size_t sep = fileName.find_last_of("\\/");
        if (sep != std::string::npos)
            fileName = fileName.substr(sep + 1, fileName.size() - sep - 1);
        size_t dot = fileName.find_last_of(".");
        if (dot != std::string::npos) {
            sampleName = fileName.substr(0, dot);
            if (fileName.substr(dot, fileName.size() - dot) == ".bam") {
                uTime.setTime();
                cerr << "Parsing sample: " << sampleName;
                fflush(NULL);
                readFactory.getSamples().push_back(sampleName);
                count = readFactory.processReadsFromBAM(bamFileName, sampleName, onlyProperlyPaired, minMAPQ, minOverlap);
                cerr << " " << count << " reads processed in " << uTime.getElapseTimeSec() << " seconds" << endl;
                fflush(NULL);
            }
        }
    }

    cerr << "Printing results" << endl;
    fflush(NULL);
    readFactory.printResults(singleFile, extendedOutput, all_feat);

    cerr << "Total time: " << uTime.getTotalTimeSec() << " seconds" << endl;
    return 0;
}