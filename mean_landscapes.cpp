#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <queue>
#include <cmath>

#include "tbb/concurrent_vector.h"

using namespace std;

#define infinity (1e30)
#define epsilon 1e-10

double get_value_at_point(const std::vector<std::pair<double, double>> &persistence_landscape_x, int position, double point) {
    if (position <= 0 || position >= persistence_landscape_x.size()) {
        return infinity;
    }
    double prev = persistence_landscape_x[position - 1].first;
    double current = persistence_landscape_x[position].first;
    if (abs(prev - current) < epsilon) {
        return (persistence_landscape_x[position - 1].second + persistence_landscape_x[position].second) / 2.0;
    }
    double l = persistence_landscape_x[position - 1].second, r = persistence_landscape_x[position].second;
    double diff_heights = r - l;
    double result = l + (point - prev) / (current - prev) * diff_heights;
    return result;
}

std::vector<std::pair<double, double>> get_mean_persistence(tbb::concurrent_vector<std::vector<std::pair<double, double>>> &persistence_landscapes_x) {
    std::vector<std::pair<double, double>> v;
    priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, greater<>> lowest_points_by_x;
    std::vector<int> pointers_by_x(persistence_landscapes_x.size(), 0);

    for (int i = 0; i < persistence_landscapes_x.size(); ++i) {
        sort(persistence_landscapes_x[i].begin(), persistence_landscapes_x[i].end());
        if (!persistence_landscapes_x[i].empty()) {
            lowest_points_by_x.push({persistence_landscapes_x[i][0].first, i});
        }
    }
    while (!lowest_points_by_x.empty()) {
        auto triple = lowest_points_by_x.top();
        double sum = 0;
        lowest_points_by_x.pop();
        sum += persistence_landscapes_x[triple.second][pointers_by_x[triple.second]].second;
        ++pointers_by_x[triple.second];
        if (pointers_by_x[triple.second] < persistence_landscapes_x[triple.second].size()) {
            lowest_points_by_x.push(
                    {persistence_landscapes_x[triple.second][pointers_by_x[triple.second]].first, triple.second});
        }
        int counter = 1;
        for (int i = 0; i < persistence_landscapes_x.size(); ++i) {
            if (i != triple.second) {
                double current_value = get_value_at_point(persistence_landscapes_x[i],
                                                          pointers_by_x[i], triple.first);
                if (current_value < infinity - 1) {
                    sum = (sum * counter + current_value) / (counter + 1);
                    ++counter;
                }
            }
        }
        if (counter > -1) {
            v.push_back({triple.first, sum});
        }
    }
    return v;
}