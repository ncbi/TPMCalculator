/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   DiffExpIR.cpp
 * Author: veraalva
 * 
 * Created on September 6, 2017, 9:50 AM
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
#include "ReadFactory.h"
#include "DiffExpIR.h"
#include "Stats.h"

using namespace std;
using namespace ngs;
using namespace genome;

void DiffExpIR::calculateDiffExpIR(ReadFactory& readFactory, std::vector<std::string> samples, std::string method, bool useFDR) {
    stats::WilcoxTest wTest;
    stats::TTest ttest;
    stats::FDRCorrection fdrCorrection;
    SPtrChromosomeNGS c;
    SPtrGeneNGS g;
    SPtrIsoformNGS i;
    SPtrFeatureNGS f;

    for (auto cIt : readFactory.getGenomeFactory().getChromosomes()) {
        c = cIt.second;
//        cout << "Chromosome: " << c->getId() << endl;
        for (auto it : c->getGenes()) {
            g = it;
//            cout << "Gene: " << g->getId();
//            fflush(NULL);
            if (g->isProcessed()) {
//                cout << "\tProcessed" << endl;
//                fflush(NULL);
                double e11_TPM, e12_TPM, e21_TPM, e22_TPM;
                int e1_count, e2_count;
                e11_TPM = e21_TPM = 0.0;
                e12_TPM = e22_TPM = 0.0;
                e1_count = e2_count = 0;
//                cout << "\tUniquefeatures: " << g->getUniquefeatures().size() << endl;
//                fflush(NULL);
                for (auto fIt = g->getUniquefeatures().begin(); fIt != g->getUniquefeatures().end(); ++fIt) {
                    f = *fIt;
//                    cout << "\t\tFeature: " << f << endl;
//                    fflush(NULL);
                    if (f->getType() == "exon") {
                        e11_TPM = e21_TPM = 0.0;
                        e1_count = e2_count = 0;
                        for (auto s : readFactory.getSamples()) {
                            try {
                                SPtrSampleData sd = f->getData().getSampleData(s);
                                if (s.compare(0, samples[0].size(), samples[0]) == 0) {
                                    e11_TPM += sd->getTPM();
                                    e1_count++;
                                } else if (s.compare(0, samples[1].size(), samples[1]) == 0) {
                                    e21_TPM += sd->getTPM();
                                    e2_count++;
                                }
                            } catch (exceptions::NotFoundException) {
                            }
                        }
                        e11_TPM = e11_TPM / static_cast<double> (e1_count);
                        e21_TPM = e21_TPM / static_cast<double> (e2_count);
                    } else if (f->getType() == "intron") {
                        vector<double> x;
                        double x_sum = 0.0;
                        vector<double> y;
                        double y_sum = 0.0;
                        for (auto s : readFactory.getSamples()) {
                            try {
                                SPtrSampleData sd = f->getData().getSampleData(s);
                                if (s.compare(0, samples[0].size(), samples[0]) == 0) {
                                    double tpm = sd->getTPM();
                                    if (tpm < 10E-5) tpm = 10E-5;
                                    x.push_back(tpm);
                                    x_sum += tpm;
                                } else if (s.compare(0, samples[1].size(), samples[1]) == 0) {
                                    double tpm = sd->getTPM();
                                    if (tpm < 10E-5) tpm = 10E-5;
                                    y.push_back(tpm);
                                    y_sum += tpm;
                                }
                            } catch (exceptions::NotFoundException) {
                            }
                        }
                        for (auto eIt = fIt; eIt != g->getUniquefeatures().end(); ++eIt) {
                            if ((*eIt)->getType() == "exon") {
                                e12_TPM = e22_TPM = 0.0;
                                e1_count = e2_count = 0;
                                for (auto s : readFactory.getSamples()) {
                                    try {
                                        SPtrSampleData sd = f->getData().getSampleData(s);
                                        if (s.compare(0, samples[0].size(), samples[0]) == 0) {
                                            e12_TPM += sd->getTPM();
                                            e1_count++;
                                        } else if (s.compare(0, samples[1].size(), samples[1]) == 0) {
                                            e22_TPM += sd->getTPM();
                                            e2_count++;
                                        }
                                    } catch (exceptions::NotFoundException) {
                                    }
                                }
                                e12_TPM = e12_TPM / static_cast<double> (e1_count);
                                e22_TPM = e22_TPM / static_cast<double> (e2_count);
                                break;
                            }
                        }
                        x_sum = x_sum / x.size();
                        y_sum = y_sum / y.size();
                        double r1 = std::log2(x_sum / (e11_TPM + e12_TPM));
                        double r2 = std::log2(y_sum / (e21_TPM + e22_TPM));
                        double p;
                        if (x.size() != 0 && y.size() != 0) {    
                            if (method == "ttest") {
                                p = ttest.pvalue(x, y);
                            } else {
                                p = wTest.pvalue(x, y);
                            }
                            if (!std::isnan(p)) {
//                                cout << "\t\t\tPValue: " << p << " R: " << r1 << " " << r2 << " Mean: " << x_sum << " " << y_sum << " log2: " << std::log2(x_sum / y_sum) << endl;
                                pvalue.push_back(p);
                                SptrDiffExpIntron d = std::make_shared<DiffExpIntron>(DiffExpIntron(make_pair(r1, r2), g, f, c->getId(), p, std::log2(x_sum / y_sum), x_sum, y_sum));
                                diffexpIRdata.push_back(d);
                            }
                        }
                    }
                }
            } else {
//                cout << "\tNo processed" << endl;
                fflush(NULL);
            }
        }
    }

    if (useFDR)
        pvalue = fdrCorrection.fdr_correction(pvalue);
}

void DiffExpIR::printDiffExpIR(std::string output_name, double fc_cutoff, double pvalue_cutoff, double r_cutoff) {
    ofstream out_file;

    out_file.open(output_name);
    out_file << "GeneId\tChr\tStart\tEnd\tIntron_Start\tIntron_End\tLog2TPMRatio\tTPM_1\tTPM_2\tminusLog10PValue\tPValue\tRValue_1\tRValue_2" << endl;
    for (unsigned int ind = 0; ind != pvalue.size(); ind++) {
        SptrDiffExpIntron d = diffexpIRdata[ind];

        if (std::isfinite(d->getRvalueFirst()) && std::isfinite(d->getRvalueSecond())) {
            if (d->getRvalueFirst() >= r_cutoff || d->getRvalueSecond() >= r_cutoff) {
                if (pvalue[ind] <= pvalue_cutoff && std::fabs(d->getLog2TPMRatio()) >= fc_cutoff) {
                    double minusLog10PValue = -1.0 * std::log10(pvalue[ind]);
                    SPtrGeneNGS g = d->getGene();
                    SPtrFeatureNGS f = d->getIntron();
                    out_file << g->getId() << "\t" << d->getChr() << "\t" << (g->getStart() + 1) << "\t" << (g->getEnd() + 1) << "\t";
                    out_file << (f->getStart() + 1) << "\t" << (f->getEnd() + 1) << "\t";
                    out_file << d->getLog2TPMRatio() << "\t" << d->getTPM_1() << "\t" << d->getTPM_2() << "\t";
                    out_file << minusLog10PValue << "\t" << pvalue[ind] << "\t" << d->getRvalueFirst() << "\t" << d->getRvalueSecond() << endl;
                }
            }
        }
    }
    out_file.close();
}

