#pragma once

#include <vector>
#include <cmath>


bool compare (const tbb::concurrent_vector<std::vector<std::pair<double, double>>>& lhs, const tbb::concurrent_vector<std::vector<std::pair<double, double>>>& rhs) {

    double special_e = 1e-10;
    double special_E = 1e30;

    if (lhs.size() != rhs.size()) {
        return false;
    }


    for (int i = 0; i < lhs.size(); ++i) {
        if (lhs[i].size() != rhs[i].size()) {
            return false;
        }
    }

    for (int i = 0; i < lhs.size(); ++i) {
        for (int j = 0; j < lhs[i].size(); ++i) {
            std::pair<double, double> l = lhs[i][j];
            std::pair<double, double> r = rhs[i][j];
            if (abs(l.first - r.first) > special_e) {
                return false;
            }
            if (abs(l.second) > special_E && abs(r.second) > special_E) {
                continue;
            }
            if ((abs(l.second) > special_E && abs(r.second) < special_E) || (abs(l.second) < special_E && abs(r.second) > special_E)) {
                return false;
            }
            if (abs(l.second - r.second) > special_e) {
                return false;
            }
        }
    }

    return true;
}