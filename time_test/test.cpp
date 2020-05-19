
int num_ = 10;

#include <atomic>
#include "..source/algorithm/landscape_a.h"
#include "..source/ripser/landscape_r.h"
#include "../util/compare.cpp"
#include <iostream>


using namespace smpl;

int main() {
    std::string path = "..dataset/magnetometer/s50.txt";
    double radii = 1e11;
    std::cout << "\n\nTest " << path << " diagram with r = 1e11" << std::endl;
    std::vector<double> v;

    tbb::concurrent_vector<std::vector<std::pair<double, double>>> diagram_a;
    tbb::concurrent_vector<std::vector<std::pair<double, double>>> diagram_r;

    for (int i = 0; i < 1; ++i) {
        double time1 = main_ripser(diagram_r, path, "/Users/leonardbee/Desktop/dataset/tore/sampled_persistence",
                                   2, radii, true, 1, 1, 1, true);
        double time2 = main_algorithm(diagram_a, path, "/Users/leonardbee/Desktop/dataset/tore/sampled_persistence",
                                      2, radii, true, 1, 1, 1, true);
    }

    auto p = get_M_D(v);
    auto res = compare(diagram_a, diagram_r);
    std::cout << (res) << std::endl;
    return 0;
}