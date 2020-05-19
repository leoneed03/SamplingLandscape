#define BOOST_TEST_MODULE testSuiteCalculator

int num_ = 3;

#include <boost/test/included/unit_test.hpp>
#include "source/algorithm/landscape_a.h"
#include <iostream>

using namespace smpl;

BOOST_AUTO_TEST_SUITE(testSuiteCalculator)
/*
    BOOST_AUTO_TEST_CASE(testCalculator0) {

            double rad = 1e11;
            std::cout << "\n\nTest sampled diagram with r = " << rad << std::endl;
            std::vector<double> v;
            tbb::concurrent_vector<std::vector<std::pair<double, double>>> diagram;
            for (int i = 0; i < 1; ++i) {
                double time = main_algorithm(diagram, "dataset/magnetometer/s50.txt", "/Users/leonardbee/Desktop/dataset/tore/sampled_persistence",
                                          2, rad, true, 1, 1, 1, true);
//                double time2 =
                v.push_back(time);
            }
            auto p = get_M_D(v);
            std::cout << "\n E = " << p.first << " D = " << p.second << '\n';

            for (int i = 0; i < diagram.size(); ++i) {
                std::cout << "\nTested " << i << std::endl;
                for (const auto& e: diagram[i]) {
                    std::cout << e.first << ' ' << e.second << std::endl;
                }

            }

//            exit(0); ////EXIT HERE

    }*/

/*

    BOOST_AUTO_TEST_CASE(testCalculator1) {

            std::cout << "\n\nTest sampled diagram with r = 8" << std::endl;
            std::vector<double> v, av;
            tbb::concurrent_vector<std::vector<std::pair<double, double>>> diagram;
            for (int i = 0; i < 1; ++i) {
                std::clog << "\n\nSTARTED " << i << std::endl;
                double time = main_algorithm(diagram, "dataset/figures/dots50_no_number.txt", "",
                                             2, 800, true, 4,8,0.4);

                std::clog << "\n\nCOMPUTED " << i << std::endl;
                v.push_back(time);
                av.push_back(matrix_size_cntr);
            }
            auto p = get_M_D(v);
            auto total_size = get_M_D(av);
            std::cout << "\n E = " << p.first << " D = " << p.second << '\n';
            std::cout << "\n size E = " << total_size.first << " size D = " << total_size.second << '\n';


    }*/



    BOOST_AUTO_TEST_CASE(testCalculator2) {
        std::cout << "\n\nTest sampled diagram with r = 0.5" << std::endl;
        std::vector<double> v, av;
        tbb::concurrent_vector<std::vector<std::pair<double, double>>> diagram;
        for (int i = 0; i < 1 /*num_*/ ; ++i) {
               std::clog << "\n\nSTARTED " << i << std::endl;
            double time = main_algorithm(diagram, "dataset/figures/human500.txt", "",
                        2, 0.5, true, 4, 8, 0.4);
            
            std::clog << "\n\nCOMPUTED " << i << std::endl;
            v.push_back(time);
            av.push_back(matrix_size_cntr);
        }
        auto p = get_M_D(v);
        auto total_size = get_M_D(av);
        std::cout << "\n E = " << p.first << " D = " << p.second << '\n';
        std::cout << "\n size E = " << total_size.first << " size D = " << total_size.second << '\n';
        
    }

BOOST_AUTO_TEST_SUITE_END()
