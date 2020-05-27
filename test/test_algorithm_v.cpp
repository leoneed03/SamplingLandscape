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

int num_test_diagram_1 = 10;


TEST(DiagramComputationFigures, bunny_basic) {

    std::string path = "../dataset/figures/bunny500.txt";
    double radii = 0.5;

    std::vector<double> v;

    tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> diagram_a;
//    tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> diagram_r;

    for (int i = 0; i < 2; ++i) {
//        double time1 = main_ripser(diagram_r, path, "", 2, radii, 4, 10, 0.3, true);
        double time2 = main_algorithm(diagram_a, path, "", 2, radii, 4, 10, 0.2, true);
    }

//    auto p = get_M_D(v);
    ASSERT_EQ(true, true);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
