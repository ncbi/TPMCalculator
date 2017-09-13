/* 
 * File:   Stats.cpp
 * Author: veraalva
 * 
 * Created on June 9, 2017, 4:22 PM
 */

/* 
 * File:   FastaFactory.cpp
 * Author: veraalva
 * 
 * Created on February 10, 2016, 3:41 PM
 */

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <memory>
#include <cmath>

#include "bmath.h"
#include "Global.h"
#include "Exceptions.h"
#include "Stats.h"

using namespace std;
using namespace stats;

vector<double> vector_rank(vector<double> r) {
    long unsigned int j = 0;
    vector<int> in;
    vector<double> rk = r;
    std::sort(rk.begin(), rk.end());
    for (auto it = rk.begin(); it != rk.end(); ++it) {
        for (auto it1 = r.begin(); it1 != r.end(); ++it1) {
            if (*it == *it1) {
                int index = it1 - r.begin();
                if (std::find(in.begin(), in.end(), index) == in.end())
                    in.push_back(it1 - r.begin());
            }
        }
    }
    for (long unsigned int i = 0; i < r.size(); i = j + 1) {
        j = i;
        while ((j < r.size() - 1) && r[in[j]] == r[in[j + 1]]) j++;
        for (long unsigned int k = i; k <= j; k++)
            rk[in[k]] = (i + j + 2) / 2.;
    }
    return rk;
}

double sum_nties(vector<double> r, vector<double> r_unique) {
    double s = 0.0;
    for (long unsigned int i = 0; i < r_unique.size(); i++) {
        double c = 0;
        for (long unsigned int j = 0; j < r.size(); j++) {
            if (r_unique[i] == r[j]) c++;
        }
        s += (std::pow(c, 3) - c);
    }
    return s;
}

double sign(double x) {
    if (x > 0) return 1;
    if (x < 0) return -1;
    return 0;

}

double Stats::variance(std::vector<double>& x, double x_mean) {
    double var = 0;
    for (auto e : x) {
        var += ((e - x_mean)*(e - x_mean));
    }
    return (var / static_cast<double> (x.size() - 1));
}

double Stats::variance(std::vector<double>& x) {
    double mx = 0.0;
    for (auto e : x) {
        mx += e;
    }
    mx = mx / static_cast<double> (x.size());
    return variance(x, mx);
}

double TTest::pvalue(std::vector<double>& x, std::vector<double>& y) {
    if (x.empty() or y.empty()) {
        throw new exceptions::EmptyDatasetException("Your dataset is empty");
    }
    bool equal = (x.size() == y.size()) ? true : false;
    double pval = NAN;
    double mx = 0.0;
    double my = 0.0;
    double df, stderr;
    Stats stats;

    for (auto e : x) {
        mx += e;
    }
    mx = mx / x.size();
    for (auto e : y) {
        my += e;
    }
    my = my / y.size();

    double vx = stats.variance(x);
    double vy = stats.variance(y);

    if (equal) {
        df = static_cast<double> (x.size() + y.size() - 2);
        double v = 0;
        if (x.size() > 1) v = v + static_cast<double> (x.size() - 1) * vx;
        if (y.size() > 1) v = v + static_cast<double> (y.size() - 1) * vy;
        v = v / df;
        stderr = std::sqrt(v * (1 / static_cast<double> (x.size()) + 1 / static_cast<double> (y.size())));
    } else {
        double stderrx = std::sqrt(vx / static_cast<double> (x.size()));
        double stderry = std::sqrt(vy / static_cast<double> (y.size()));
        stderr = std::sqrt(std::pow(stderrx, 2) + std::pow(stderry, 2));
        df = std::pow(stderr, 4) / (std::pow(stderrx, 4) / static_cast<double> (x.size() - 1) + std::pow(stderry, 4) / static_cast<double> (y.size() - 1));
    }
    double tstat = (mx - my) / stderr;
    pval = 2 * pt(-std::abs(tstat), df, 1, 0);
    return pval;
}

double WilcoxTest::pvalue(std::vector<double>& x, std::vector<double>& y) {
    if (x.empty() or y.empty()) {
        throw new exceptions::EmptyDatasetException("Your dataset is empty");
    }
    bool exact = (x.size() < 50) && (y.size() < 50);
    double p = NAN;
    double stats = 0.0;

    vector<double> r(x.begin(), x.end());
    for (auto it = y.begin(); it != y.end(); ++it) {
        r.push_back(*it);
    }
    r = vector_rank(r);
    for (long unsigned int i = 0; i < x.size(); i++) stats += r[i];
    stats -= static_cast<double> (x.size()) * (static_cast<double> (x.size()) + 1) / 2;

    vector<double> r_unique(r.begin(), r.end());
    std::sort(r_unique.begin(), r_unique.end());
    auto last = std::unique(r_unique.begin(), r_unique.end());
    r_unique.erase(last, r_unique.end());
    bool ties = (r.size() != r_unique.size());

    //    cout << "Len(r): " << r.size() << endl;
    //    cout << "Len(unique(r)):: " << r_unique.size() << endl;
    //    cout << "Is exact:: " << exact << endl;
    //    cout << "TIES: " << ties << endl;

    if (exact && !ties) {
        //        printf("exact && !TIES\n");
        if (stats > (static_cast<double> (x.size()) * static_cast<double> (y.size()) / 2.0)) {
            p = pwilcox(stats - 1, static_cast<double> (x.size()), static_cast<double> (y.size()), false, false);
        } else {
            p = pwilcox(stats, static_cast<double> (x.size()), static_cast<double> (y.size()), true, false);
        }
        if (std::isnan(p)) return NAN;
        p = min(2 * p, 1);
    } else {
        //        printf("!exact && !TIES\n");
        double z = stats - static_cast<double> (x.size()) * static_cast<double> (y.size()) / 2.0;
        double sigma = sqrt((static_cast<double> (x.size()) * static_cast<double> (y.size()) / 12) \
            * ((static_cast<double> (x.size()) + static_cast<double> (y.size()) + 1) - \
            sum_nties(r, r_unique) / ((static_cast<double> (x.size()) + static_cast<double> (y.size())) \
            * (static_cast<double> (x.size()) + static_cast<double> (y.size()) - 1))));
        double correction = sign(z) * 0.5;
        //        printf("Z: %f\n", z);
        //        printf("SIGMA: %f\n", sigma);
        //        printf("CORRECTION: %f\n", correction);
        z = (z - correction) / sigma;
        //        printf("Z: %f\n", z);
        if (std::isnan(z)) return NAN;
        p = 2 * min(pnorm5(z, 0, 1, true, false), pnorm5(z, 0, 1, false, false));
    }
    return p;
}

vector<double> cummin(vector<double> & c) {
    vector<double> cmin;
    double min = INFINITY;
    for (auto it = c.begin(); it != c.end(); ++it) {
        if (std::isnan(*it) || std::isnan(min))
            min = min + *it; /* propagate NA and NaN */
        else
            min = (min < *it) ? min : *it;
        cmin.push_back(min);
    }
    return cmin;
}

std::vector<int> sorted_order(std::vector<double> &x, bool decreasing) {
    std::vector<int> y(x.size());
    std::size_t n(0);
    std::generate(std::begin(y), std::end(y), [&] {
        return n++; });

    if (decreasing) {
        std::sort(std::begin(y),
                std::end(y),
                [&](int i1, int i2) {
                    return x[i1] > x[i2]; });
    } else {
        std::sort(std::begin(y),
                std::end(y),
                [&](int i1, int i2) {
                    return x[i1] < x[i2]; });
    }
    return y;
}

std::vector<int> sorted_order(std::vector<int> &x, bool decreasing) {
    std::vector<int> y(x.size());
    std::size_t n(0);
    std::generate(std::begin(y), std::end(y), [&] {
        return n++; });

    if (decreasing) {
        std::sort(std::begin(y),
                std::end(y),
                [&](int i1, int i2) {
                    return x[i1] > x[i2]; });
    } else {
        std::sort(std::begin(y),
                std::end(y),
                [&](int i1, int i2) {
                    return x[i1] < x[i2]; });
    }
    return y;
}

/**
 * Calculates FDR correction.
 * It assumes no NAN are included 
 * 
 * @param c Vector of floats (P-Values)
 * @return 
 */
vector<double> FDRCorrection::fdr_correction(vector<double> & c) {
    vector<double> fdr;
    vector<double> f;
    int n = c.size();
    vector<int> c_indexes_sorted = sorted_order(c, true);
    vector<int> rc_indexes_sorted = sorted_order(c_indexes_sorted, false);

    int i = c.size();
    for (auto v : c_indexes_sorted) {
        double m = static_cast<double> (n) / static_cast<double> (i) * c[v];
        m = (m > 1) ? 1 : m;
        f.push_back(m);
        i--;
    }
    for (auto v : rc_indexes_sorted) {
        fdr.push_back(f[v]);
    }
    return fdr;
}


