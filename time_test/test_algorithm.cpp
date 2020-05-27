int num_ = 3;

//#include <boost/test/included/unit_test.hpp>
#include <gtest/gtest.h>
//#include "/Users/leonardbee/CLionProjects/SubsamplingMethodsForPersistenceLandscape1/ripser.cpp"
//#include "/Users/leonardbee/CLionProjects/SubsamplingMethodsForPersistenceLandscape1/mean_landscapes.cpp"
//#include "../include/gudhi/landscape_g.h"
#include "../include/mean_landscapes/mean_landscapes.h"
#include "../include/algorithm/landscape_a.h"
#include <iostream>



using namespace smpl;


//TEST(DiagramComputationFigures, bunny_basic) {
//
//    std::string path = "../dataset/figures/bunny500.txt";
//    double radii = 0.2;
//
//    std::vector<double> v;
//
//    tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> diagram_a;
//    tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> diagram_r;
//
//    for (int i = 0; i < 1; ++i) {
//        double time1 = main_ripser(diagram_r, path, "", 2, radii, 1, 1, 1, true);
//        double time2 = main_algorithm(diagram_a, path, "", 2, radii, 1, 1, 1, true);
//    }
//
//    auto p = get_M_D(v);
//    auto res = compare(diagram_a[0], diagram_r[0]);
//    ASSERT_EQ(res, true);
//}

TEST(LandscapeComputationAlgorithm, random_500) {
    std::cout << "\n\nTest sampled diagram with r =  diam" << std::endl;
    std::vector<double> v;
    tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> apd;
    for (int i = 0; i < num_; ++i) {
        double time = landscape_algorithm("../dataset/random_500.txt", "results", 2, 0.5, 4, 10, 0.4);
        v.push_back(time);
    }
    auto p = get_M_D(v);
    std::cout << "\n E = " << p.first << " D = " << p.second << '\n';
}
TEST(LandscapeComputationAlgorithm, magnetometer_500) {
    std::cout << "\n\nTest sampled diagram with r =  diam" << std::endl;
    std::vector<double> v;
    tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> apd;
    for (int i = 0; i < num_; ++i) {
        double time = landscape_algorithm("../dataset/magnetometer/basketball/p2_1000.txt", "results", 2, 0.5, 4, 10, 0.4);
        v.push_back(time);
    }
    auto p = get_M_D(v);
    std::cout << "\n E = " << p.first << " D = " << p.second << '\n';
}

/*
TEST(LandscapeComputationAlgorithm, human_500) {
    std::cout << "\n\nTest sampled diagram with r =  diam" << std::endl;
    std::vector<double> v;
    tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> apd;
    for (int i = 0; i < num_; ++i) {
        double time = landscape_algorithm("../dataset/figures/human500.txt", "results", 2, 0.5, 4, 10, 0.4);
        v.push_back(time);
    }
    auto p = get_M_D(v);
    std::cout << "\n E = " << p.first << " D = " << p.second << '\n';
}
*/

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
