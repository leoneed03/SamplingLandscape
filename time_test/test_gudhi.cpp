#define BOOST_TEST_MODULE testSuiteCalculator

int num_ = 3;

//#include <boost/test/included/unit_test.hpp>
#include <gtest/gtest.h>
//#include "/Users/leonardbee/CLionProjects/SubsamplingMethodsForPersistenceLandscape1/ripser.cpp"
//#include "/Users/leonardbee/CLionProjects/SubsamplingMethodsForPersistenceLandscape1/mean_landscapes.cpp"
//#include "../include/gudhi/landscape_g.h"
#include "../include/mean_landscapes/mean_landscapes.h"
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

TEST(LandscapeComputationGudhi, random_500) {
    std::cout << "\n\nTest sampled diagram with r =  diam" << std::endl;
    std::vector<double> v;
    tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> apd;
    for (int i = 0; i < num_; ++i) {
        double time = landscape_with_gudhi("../dataset/random_500.txt", "results", 2, 0.5, 4, 10, 0.4);
        v.push_back(time);
    }
    auto p = get_M_D(v);
    std::cout << "\n E = " << p.first << " D = " << p.second << '\n';
}
TEST(LandscapeComputationGudhi, magnetometer_500) {
    std::cout << "\n\nTest sampled diagram with r =  diam" << std::endl;
    std::vector<double> v;
    tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> apd;
    for (int i = 0; i < num_; ++i) {
        double time = landscape_with_gudhi("../dataset/magnetometer/basketball/p2_1000.txt", "results", 2, 0.5, 4, 10, 0.4);
        v.push_back(time);
    }
    auto p = get_M_D(v);
    std::cout << "\n E = " << p.first << " D = " << p.second << '\n';
}

TEST(LandscapeComputationGudhi, human_500) {
    std::cout << "\n\nTest sampled diagram with r =  diam" << std::endl;
    std::vector<double> v;
    tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> apd;
    for (int i = 0; i < num_; ++i) {
        double time = landscape_with_gudhi("../dataset/figures/human500.txt", "results", 2, 0.5, 4, 10, 0.4);
        v.push_back(time);
    }

//    std::cout << "S_amples " << apd.size() << std::endl;
//    for (const auto& e: apd) {
//        std::cout << "            dims " << e.size() << std::endl;
//        for (const auto& a: e) {
//            std::cout << "                         intervals " << a.size() << std::endl;
//        }
//    }
    auto p = get_M_D(v);
    std::cout << "\n E = " << p.first << " D = " << p.second << '\n';
}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

/*
    BOOST_AUTO_TEST_CASE(testCalculator4) {
        std::ofstream out;
        out.open("loggudhi.txt", std::ios::app);
        out << "\n\n\nnew loggudhi.n r=diam full ";
        std::cout << "\n\nTest full diagram with r = 0.5" << std::endl;
        std::vector<double> v;
        for (int i = 0; i < num_; ++i) {
            double time = landscape_gudhi("../dataset/figures/human500.txt", "results",
                        2, 1e11, 1, 1, 1);
            v.push_back(time);
            out << time << ' ';
        }
        auto p = get_M_D(v);
        out << "\n E = " << p.first << " D = " << p.second << '\n';
        std::cout << "\n E = " << p.first << " D = " << p.second << '\n';
        out.close();
//        BOOST_CHECK_EQUAL(s(6,6), 12);
    }
 */
