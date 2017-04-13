#include <iostream>
#include <ctime>
#include <cstdint>
#include <set>
#include <unordered_map>
#include <algorithm>
#include <memory>

#include "api/BamReader.h"

#include "Global.h"
#include "Exceptions.h"
#include "TimeUtils.h"
#include "bstring.h"
#include "TextParser.h"
#include "GenomeFactory.h"
#include "ReadFactory.h"

using namespace std;
using namespace ngs;
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
    cerr << "-c    Control sample identification prefix. This is a prefix pattern in the control BAM file names\n";
    cerr << "-t    Treated sample identification prefix. This is a prefix pattern in the treated BAM file names\n";


    cerr << "\n********************************************************************************\n";
    cerr << "\n                        Roberto Vera Alvarez, PhD\n";
    cerr << "            Emails: veraalva@ncbi.nlm.nih.gov, r78v10a07@gmail.com\n\n";
    cerr << "********************************************************************************\n";
    exit(exit_code);
}

int main(int argc, char *argv[]) {
    int count = 0;
    TimeUtils uTime;
    string gtfFileName;
    string bamDirName;
    set<string>features = {"exon"};
    unordered_map<string, string> featuresToCreate = {
        {"exon", "intron"}
    };
    string ctrl, treated;
    ReadFactory readFactory;
    
    if (argc == 1){
        print_usage(argv[0], 0);
    }

    for (int i = 1; i < argc; i++) {
        string option(argv[i]);
        if (option.compare(0, 1, "-") == 0 && option.compare(1, 1, "-") != 0 && option.size() == 2) {
            if (option.compare(1, 1, "h") == 0) {
                print_usage(argv[0], 0);
            } else if (option.compare(1, 1, "v") == 0) {
                Global::instance()->setVerbose(1);
            } else if (option.compare(1, 1, "d") == 0) {
                Global::instance()->setVerbose(3);
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
            } else if (option.compare(1, 1, "c") == 0) {
                i++;
                if (i < argc) {
                    ctrl = argv[i];
                    if (ctrl.compare(0, 1, "-") == 0) {
                        cerr << "Option c require an argument" << endl;
                        print_usage(argv[0], -1);
                    }

                } else {
                    cerr << "Option c require an argument" << endl;
                    print_usage(argv[0], -1);
                }
            } else if (option.compare(1, 1, "t") == 0) {
                i++;
                if (i < argc) {
                    treated = argv[i];
                    if (treated.compare(0, 1, "-") == 0) {
                        cerr << "Option t require an argument" << endl;
                        print_usage(argv[0], -1);
                    }

                } else {
                    cerr << "Option t require an argument" << endl;
                    print_usage(argv[0], -1);
                }
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
        cerr << "\nGTF file is required. See -g option" << endl;
        print_usage(argv[0], -1);
    }

    if (bamDirName.empty()) {
        cerr << "\nDirectory with the BAM files is required. See -d option" << endl;
        print_usage(argv[0], -1);
    }

    if (ctrl.empty()) {
        cerr << "\nControl identification string is required. See -c option" << endl;
        print_usage(argv[0], -1);
    }

    if (treated.empty()) {
        cerr << "\nTreated identification string is required. See -d option" << endl;
        print_usage(argv[0], -1);
    }
    vector<string> groupPrefix = {ctrl, treated};


    uTime.setTime();
    cerr << "Reading GTF file ... " << endl;
    readFactory.getGenomeFactory().processGTFFile(gtfFileName, "gene_name", "transcript_id", features, featuresToCreate);
    cerr << "Done in " << uTime.getElapseTimeSec() << " seconds" << endl;

    uTime.setTime();
    cerr << "Parsing BAM file" << endl;
    count = readFactory.processBAMSAMFromDir(bamDirName, groupPrefix);
    cerr << count << " reads processed in " << uTime.getElapseTimeSec() << " seconds" << endl;

    cerr << "Printing results" << endl;
    readFactory.printResults(groupPrefix);

    cerr << "Total time: " << uTime.getTotalTimeSec() << " seconds" << endl;
    return 0;
}