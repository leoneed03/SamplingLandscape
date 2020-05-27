#pragma once

#include "tbb/concurrent_vector.h"
#include <vector>
#include <iostream>
#include <cmath>

#define DEBUG_COMPARER_FLAG false

bool compare(tbb::concurrent_vector<std::vector<std::pair<double, double>>> &lhs,
             tbb::concurrent_vector<std::vector<std::pair<double, double>>> &rhs) {

    double special_e = 1e-10;
    double special_E = 1e30;

    if (DEBUG_COMPARER_FLAG) {
        std::cout << "here 0" << std::endl;
    }
    while (lhs.size() < rhs.size()) {
        lhs.push_back(std::vector<std::pair<double, double>>(0));
    }
    while (lhs.size() > rhs.size()) {
        rhs.push_back(std::vector<std::pair<double, double>>(0));
    }

    if (DEBUG_COMPARER_FLAG) {
        std::cout << "here 1" << std::endl;
    }

    for (int i = 0; i < lhs.size(); ++i) {
        if (lhs[i].size() != rhs[i].size()) {
            return false;
        }
    }

    if (DEBUG_COMPARER_FLAG) {
        std::cout << "here 2" << std::endl;
    }


    for (int i = 0; i < lhs.size(); ++i) {

        if (DEBUG_COMPARER_FLAG) {
            std::cout << "here first " << i << std::endl;
        }
        for (int j = 0; j < lhs[i].size(); ++j) {

            if (DEBUG_COMPARER_FLAG) {
                std::cout << "                       here second " << j << std::endl;
            }
            std::pair<double, double> l = lhs[i][j];
            std::pair<double, double> r = rhs[i][j];
            if (abs(l.first - r.first) > special_e) {
                return false;
            }
            if (abs(l.second) > special_E && abs(r.second) > special_E) {
                continue;
            }
            if ((abs(l.second) > special_E && abs(r.second) < special_E) ||
                (abs(l.second) < special_E && abs(r.second) > special_E)) {
                return false;
            }
            if (abs(l.second - r.second) > special_e) {
                return false;
            }
        }
    }

    return true;
}