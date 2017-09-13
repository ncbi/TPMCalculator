/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Stats.h
 * Author: veraalva
 *
 * Created on June 9, 2017, 4:22 PM
 */

#ifndef STATS_H
#define STATS_H

namespace stats {

    class WilcoxTest {
    public:

        WilcoxTest() {
        }

        virtual ~WilcoxTest() {
        }

        double pvalue(std::vector<double> &x, std::vector<double> &y);
    private:

    };

    class TTest {
    public:

        TTest() {
        }

        virtual ~TTest() {
        }

        double pvalue(std::vector<double> &x, std::vector<double> &y);
    private:

    };

    class Stats {
    public:

        Stats(){
        }

        virtual ~Stats() {
        }

        double variance(std::vector<double> &x);
        double variance(std::vector<double> &x, double x_mean);
    private:
    };

    class FDRCorrection {
    public:

        FDRCorrection() {
        }

        virtual ~FDRCorrection() {
        }

        std::vector<double> fdr_correction(std::vector<double> & c);
    private:
    };
}

#endif /* STATS_H */

