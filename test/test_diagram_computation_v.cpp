#include <gtest/gtest.h>


#include <atomic>
#include "../include/algorithm/landscape_a.h"
#include "../include/ripser/landscape_r.h"
#include "../include/mean_landscapes/mean_landscapes.h"
//#include <boost/test/included/unit_test.hpp>
#include "../util/compare.cpp"
#include <iostream>

#define PRINT_PAIR_NUMBER false
#define GTEST_HAS_PARAM_TEST 1

using namespace smpl;

int num_test_diagram_1 = 10;


TEST(DiagramComputationFigures, bunny_basic) {

    std::string path = "../dataset/figures/bunny500.txt";
    double radii = 0.2;

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

TEST(DiagramComputationMagnetometer, magnetometer_s2) {

    std::string path = "../dataset/magnetometer/s2.txt";
    double radii = 0.1;

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


TEST(DiagramComputationMagnetometer, magnetometer_s3) {

    std::string path = "../dataset/magnetometer/s3.txt";
    double radii = 0.1;

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


TEST(DiagramComputationMagnetometer, magnetometer_s0) {

    std::string path = "../dataset/magnetometer/s0.txt";
    double radii = 1e11;

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


TEST(DiagramComputationMagnetometer, magnetometer_s50) {

    std::string path = "../dataset/magnetometer/s50.txt";
    double radii = 1e11;

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

TEST(DiagramComputationFigures, figures_dots50) {
    for (int j = 0; j < 10; ++j) {
    
        std::cout << "ATTEMPT " << j << std::endl;
        std::string path = "../dataset/figures/dots50.txt";
        double radii = 1e11;

        std::vector<double> v;

        tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> diagram_a;
        tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> diagram_r;

        for (int i = 0; i < 1; ++i) {

            double time1 = main_ripser(diagram_r, path, "", 2, radii, 1, 1, 1, true);
            double time2 = main_algorithm(diagram_a, path, "", 2, radii, 1, 1, 1, true);
        }
        auto res = compare(diagram_a[0], diagram_r[0]);
        ASSERT_EQ(res, true);
    }
}


TEST(DiagramComputationFigures, figures_dots50_without_point) {

    std::string path = "../dataset/figures/dots50_no_number.txt";
    double radii = 1e11;

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

TEST(DiagramComputationMagnetometer, magnetometer_s50_r5) {

    std::string path = "../dataset/magnetometer/s50.txt";
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

TEST(DiagramComputationFigures, figures_dots50_r_0_5) {
    std::string path = "../dataset/figures/dots50.txt";
    double radii = 0.5;

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


TEST(DiagramComputationMagnetometer, magnetometer_s1000) {

    std::string path = "../dataset/magnetometer/s1000.txt";
    double radii = 0.5;

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


TEST(DiagramComputationFigures, figures_bunny_500_r_0_2) {

    std::string path = "../dataset/figures/bunny500.txt";
    double radii = 0.2;

    std::vector<double> v;

    tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> diagram_a;
    tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> diagram_r;

    for (int i = 0; i < 1; ++i) {
        double time1 = main_ripser(diagram_r, path, "", 2, radii, 1, 1, 1, true);
        double time2 = main_algorithm(diagram_a, path, "", 2, radii, 1, 1, 1, true);
    }

    auto p = get_M_D(v);
    bool res = true;

    res = compare(diagram_a[0], diagram_r[0]);

    ASSERT_EQ(res, true);
}

TEST(DiagramComputationFigures, figures_sphere_500_r_0_2) {

    std::string path = "../dataset/figures/sphere500.txt";
    double radii = 0.2;

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

TEST(DiagramComputationFigures, figures_tore_500_r_0_2) {

    std::string path = "../dataset/figures/tore500.txt";
    double radii = 0.2;


    std::vector<double> v;

    tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> diagram_a;
    tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> diagram_r;

    for (int i = 0; i < 1; ++i) {
        double time1 = main_ripser(diagram_r, path, "", 2, radii, 1, 1, 1, true);
        double time2 = main_algorithm(diagram_a, path, "", 2, radii, 1, 1, 1, true);
    }

    auto p = get_M_D(v);
    auto res = compare(diagram_a[0], diagram_r[0]);

    ASSERT_EQ(res, false);
}


TEST(DiagramComputationFigures, figures_human_500) {

    std::string path = "../dataset/figures/human500.txt";
    double radii = 0.2;

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

    if (PRINT_PAIR_NUMBER) {
        std::cout << "S_amples a" << diagram_a.size() << std::endl;
        for (const auto &e: diagram_a) {
            std::cout << "            dims " << e.size() << std::endl;
            for (const auto &a: e) {
                std::cout << "                         intervals " << a.size() << std::endl;
            }
        }

        std::cout << "S_amples b " << diagram_r.size() << std::endl;
        for (const auto &e: diagram_r) {
            std::cout << "            dims " << e.size() << std::endl;
            for (const auto &a: e) {
                std::cout << "                         intervals " << a.size() << std::endl;
            }
        }
    }

}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
