#define BOOST_TEST_MODULE testSuiteCalculator

int num_test_diagram_1 = 10;

#include <atomic>
#include "../include/algorithm/landscape_a.h"
#include "../include/ripser/landscape_r.h"
#include <boost/test/included/unit_test.hpp>
#include "../util/compare.cpp"
#include <iostream>

using namespace smpl;

BOOST_AUTO_TEST_SUITE(testSuiteCalculator)

    BOOST_AUTO_TEST_CASE(testDiagram5) {

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

        BOOST_CHECK_EQUAL(res, true);
    }

/*


    BOOST_AUTO_TEST_CASE(testOne2) {

        std::string path = "../dataset/magnetometer/s2.txt";
        double radii =  0.1;

        std::vector<double> v;

        tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> diagram_a;
        tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> diagram_r;

        for (int i = 0; i < 1; ++i) {
            double time1 = main_ripser(diagram_r, path, "", 2, radii, 1, 1, 1, true);
            double time2 = main_algorithm(diagram_a, path, "", 2, radii, 1, 1, 1, true);
        }

        auto p = get_M_D(v);
        auto res = compare(diagram_a, diagram_r);

        BOOST_CHECK_EQUAL(res, true);

    }

    BOOST_AUTO_TEST_CASE(testOne3) {

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

        auto res = compare(diagram_a, diagram_r);

        BOOST_CHECK_EQUAL(res, true);
    }
    

    BOOST_AUTO_TEST_CASE(testEmpty0) {

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
        auto res = compare(diagram_a, diagram_r);

        BOOST_CHECK_EQUAL(res, true);
    }



   
    BOOST_AUTO_TEST_CASE(testDiagram0) {

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
        auto res = compare(diagram_a, diagram_r);

        BOOST_CHECK_EQUAL(res, true);
    }

    BOOST_AUTO_TEST_CASE(testDiagram1) {

        std::string path = "../dataset/figures/dots50.txt";
        double radii = 1e11;
         
        std::vector<double> v;

        tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> diagram_a;
        tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> diagram_r;

        for (int i = 0; i < 1; ++i) {

            double time1 = main_ripser(diagram_r, path, "", 2, radii, 1, 1, 1, true);
            double time2 = main_algorithm(diagram_a, path, "", 2, radii, 1, 1, 1, true);
        }

        auto p = get_M_D(v);
        auto res = compare(diagram_a, diagram_r);

        BOOST_CHECK_EQUAL(res, true);
    }

    BOOST_AUTO_TEST_CASE(testDiagram2) {

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
        auto res = compare(diagram_a, diagram_r);

        BOOST_CHECK_EQUAL(res, true);
    }

    BOOST_AUTO_TEST_CASE(testDiagram3) {

        std::string path = "../dataset/magnetometer/s50.txt";
        double radii = 5;
         
        std::vector<double> v;

        tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> diagram_a;
        tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> diagram_r;

        for (int i = 0; i < 1; ++i) {
            double time1 = main_ripser(diagram_r, path, "",
                    2, radii, 1, 1, 1, true);
            double time2 = main_algorithm(diagram_a, path, "",
                    2, radii, 1, 1, 1, true);
        }

        auto p = get_M_D(v);
        auto res = compare(diagram_a, diagram_r);

        BOOST_CHECK_EQUAL(res, true);
    }

    BOOST_AUTO_TEST_CASE(testDiagram4) {

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
        auto res = compare(diagram_a, diagram_r);

        BOOST_CHECK_EQUAL(res, true);
    }


*/
   

/*
    BOOST_AUTO_TEST_CASE(testDiagram6) {

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
        auto res = compare(diagram_a, diagram_r);

        BOOST_CHECK_EQUAL(res, true);
    }

    BOOST_AUTO_TEST_CASE(testDiagram7) {

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
        auto res = compare(diagram_a, diagram_r);

        BOOST_CHECK_EQUAL(res, true);
    }

    BOOST_AUTO_TEST_CASE(testDiagram8) {

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
        auto res = compare(diagram_a, diagram_r);

        BOOST_CHECK_EQUAL(res, true);
    }


    BOOST_AUTO_TEST_CASE(testDiagram9) {

        std::string path = "../dataset/figures/human500.txt";
        double radii = 0.2;

        std::vector<double> v;

        tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> diagram_a;
        tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> diagram_r;

        for (int i = 0; i < 5; ++i) {
            double time1 = main_ripser(diagram_r, path, "", 2, radii, 1, 1, 1, true);
            double time2 = main_algorithm(diagram_a, path, "", 2, radii, 1, 1, 1, true);
        }

        auto p = get_M_D(v);
        auto res = compare(diagram_a, diagram_r);

        BOOST_CHECK_EQUAL(res, true);
    }*/


BOOST_AUTO_TEST_SUITE_END()
