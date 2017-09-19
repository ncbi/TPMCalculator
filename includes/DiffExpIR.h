/* 
 * File:   DiffExpIR.h
 * Author: veraalva
 *
 * Created on September 6, 2017, 9:50 AM
 */

#ifndef DIFFEXPIR_H
#define DIFFEXPIR_H

#include "ReadFactory.h"

namespace ngs {

    class DiffExpIntron {
    public:

        DiffExpIntron(std::pair<double, double> rvalue, SPtrGeneNGS g, SPtrFeatureNGS i, std::string chr, double pvalue, double log2TPMRatio, double TPM_1, double TPM_2) :
        rvalue(rvalue), g(g), i(i), chr(chr), pvalue(pvalue), log2TPMRatio(log2TPMRatio), TPM_1(TPM_1), TPM_2(TPM_2) {
        }

        virtual ~DiffExpIntron() {

        }

        std::string getChr() const {
            return chr;
        }

        SPtrGeneNGS& getGene() {
            return g;
        }

        SPtrFeatureNGS& getIntron() {
            return i;
        }

        double getLog2TPMRatio() {
            return log2TPMRatio;
        }

        double getPvalue() {
            return pvalue;
        }

        double getRvalueFirst() {
            return rvalue.first;
        }

        double getRvalueSecond() {
            return rvalue.second;
        }

        double getTPM_1() const {
            return TPM_1;
        }

        double getTPM_2() const {
            return TPM_2;
        }

    private:
        std::pair<double, double> rvalue;
        SPtrGeneNGS g;
        SPtrFeatureNGS i;
        std::string chr;
        double pvalue;
        double log2TPMRatio;
        double TPM_1;
        double TPM_2;
    };

    typedef std::shared_ptr<DiffExpIntron> SptrDiffExpIntron;

    class DiffExpIR {
    public:

        DiffExpIR() {
        }

        virtual ~DiffExpIR() {
        }

        void calculateDiffExpIR(ReadFactory& readFactory, std::vector<std::string> samples, std::string method, bool useFDR);

        void printDiffExpIR(std::string output_name, double fc_cutoff, double pvalue_cutoff, double r_cutoff);

    private:
        std::vector<double> pvalue;
        std::vector<SptrDiffExpIntron> diffexpIRdata;
    };

}

#endif /* DIFFEXPIR_H */

