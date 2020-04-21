#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <queue>

using namespace std;
#define infty (1e6)

double get_value_at_point_test(vector<double>& persistence_landscape_x, const vector<double>& persistence_landscape_y, double point) {
    vector<double>::iterator it = upper_bound(persistence_landscape_x.begin(), persistence_landscape_x.end(), point);
    if (it == persistence_landscape_x.begin() || it == persistence_landscape_x.end()) {
        return infty;
    }
    int position = it - persistence_landscape_x.begin();
    double prev = persistence_landscape_x[position - 1];
    double l =  persistence_landscape_y[position - 1], r = persistence_landscape_y[position];
    double diff_heights = r - l;
    double result = l + (point - prev) / (*it - prev) * diff_heights;
    return result;
}

double get_value_at_point_test(const vector<double>& persistence_landscape_x, const vector<double>& persistence_landscape_y, int position, double point) {
    if (position <= 0 || position >= persistence_landscape_x.size()) {
        return infty;
    }
    double prev = persistence_landscape_x[position - 1], current = persistence_landscape_x[position];
    double l =  persistence_landscape_y[position - 1], r = persistence_landscape_y[position];
    double diff_heights = r - l;
    double result = l + (point - prev) / (current - prev) * diff_heights;
    return result;
}
vector<pair<double, double>> get_mean_persistence_test(const vector<vector<double>>& persistence_landscapes_x, const vector<vector<double>>& persistence_landscapes_y) {
    vector<pair<double, double>> v;
    priority_queue<pair<double, int>, vector<pair<double, int>>, greater<>> lowest_points_by_x;
    vector<int> pointers_by_x (persistence_landscapes_x.size(), 0);

    for (int i = 0; i < persistence_landscapes_x.size(); ++i) {
        if (!persistence_landscapes_x[i].empty()) {
            lowest_points_by_x.push({persistence_landscapes_x[i][0], i});
        }
    }
    while (!lowest_points_by_x.empty()) {
        auto triple = lowest_points_by_x.top();
        double sum = 0;
        lowest_points_by_x.pop();
        sum += persistence_landscapes_y[triple.second][pointers_by_x[triple.second]];
        ++pointers_by_x[triple.second];
        if (pointers_by_x[triple.second] < persistence_landscapes_x[triple.second].size()) {
            lowest_points_by_x.push({persistence_landscapes_x[triple.second][pointers_by_x[triple.second]], triple.second});
        }
        int counter = 1;
        for (int i = 0; i < persistence_landscapes_x.size(); ++i) {
            if (i != triple.second) {
                double current_value = get_value_at_point_test(persistence_landscapes_x[i], persistence_landscapes_y[i],
                                                          pointers_by_x[i], triple.first);
                if (current_value < infty - 1) {
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
int test_test() {
    vector<double> x = {1, 2, 3, 4, 5, 6};
    vector<double> y = {101, 102, 103, 113, 105, 106};
    vector<double> p = {0};


//    vector<vector<double>> v_x = {{0, 4, 6, 7, 10}, {-1, 2, 5, 8}};
//    vector<vector<double>> v_y = {{4, 2, 10, 1, -5}, {20, 20, 1, 5}};
    vector<vector<double>> v_x = {{2, 5, 8, 12, 16}, {0, 5, 14, 15}};
    vector<vector<double>> v_y = {{11, 5, 7, 7, -1}, {0, 9, 13, 17}};
    int c = 0;
    cout << get_value_at_point_test(x, y, 3, 3.2) << endl;
    for (const auto& e: get_mean_persistence_test(v_x, v_y)) {
        cout << e.first << " -> " << e.second << endl;
    }
    return 0;
}


/*
 0 -> 0
2 -> 7
5 -> 6.25
6 -> 7.33333
8 -> 8.5
12 -> 9.5
14 -> 8
15 -> 1
 */