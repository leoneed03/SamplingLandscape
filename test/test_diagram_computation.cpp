#include <gtest/gtest.h>


#include <atomic>
#include "../include/algorithm/landscape_a.h"
#include "../include/ripser/landscape_r.h"
#include "../include/mean_landscapes/mean_landscapes.h"
//#include <boost/test/included/unit_test.hpp>
#include "../util/compare.cpp"
#include <iostream>

#define PRINT_PAIR_NUMBER false

using namespace smpl;


TEST(DiagramComputationMagnetometer, magnetometer_basic) {

    std::string path = "../dataset/magnetometer/s1000.txt";
    double radii = 5;

    std::vector<double> v;

    tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> diagram_a;
    tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> diagram_r;

    for (int i = 0; i < 1; ++i) {
        double time1 = main_ripser(diagram_r, path, "", 2, radii, 1, 1, 1, true);
        double time2 = main_algorithm(diagram_a, path, "", 2, radii, 1, 1, 1, true);
    }

    auto p = get_M_D(v);
    auto res = compare(diagram_a[0], diagram_r[0]);
    ASSERT_EQ(res, true);
}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
