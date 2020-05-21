#define BOOST_TEST_MODULE testSuiteCalculator

int num_ = 3;

#include <boost/test/included/unit_test.hpp>
//#include "/Users/leonardbee/CLionProjects/SubsamplingMethodsForPersistenceLandscape1/ripser.cpp"
//#include "/Users/leonardbee/CLionProjects/SubsamplingMethodsForPersistenceLandscape1/mean_landscapes.cpp"
#include "../include/gudhi/landscape_g.h"
#include <iostream>


using namespace smpl;

BOOST_AUTO_TEST_SUITE(testSuiteCalculator)
    
    BOOST_AUTO_TEST_CASE(testCalculator1) {
        std::cout << "\n\nTest sampled diagram with r = 0.5" << std::endl;
        std::vector<double> v;
        for (int i = 0; i < 1; ++i) {
            double time = landscape_gudhi("../dataset/figures/dots50_no_number.txt", "results",
                        2, 8, 1, 1, 1);
            v.push_back(time);
        }
        auto p = get_M_D(v);
        std::cout << "\n E = " << p.first << " D = " << p.second << '\n';
        
    }

    BOOST_AUTO_TEST_CASE(testCalculator2) {
        std::cout << "\n\nTest full diagram with r = 0.5" << std::endl;
        std::vector<double> v;
        for (int i = 0; i < 1; ++i) {
            double time = landscape_gudhi("../dataset/figures/human500.txt", "results",
                        2, 0.5, 1, 1, 0.3);
            v.push_back(time);
        }
        auto p = get_M_D(v);
        std::cout << "\n E = " << p.first << " D = " << p.second << '\n';
    }

    
    BOOST_AUTO_TEST_CASE(testCalculator3) {
        std::cout << "\n\nTest sampled diagram with r =  diam" << std::endl;
        std::vector<double> v;
        tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> apd;
        for (int i = 0; i < num_; ++i) {
            double time = landscape_gudhi_with_diagrams("../dataset/figures/human500.txt", "results", apd,
                        2, 0.5, 4, 10, 0.4);
            v.push_back(time);
        }

//        std::cout << "S_amples " << apd.size() << std::endl;
//        for (const auto& e: apd) {
//            std::cout << "            dims " << e.size() << std::endl;
//            for (const auto& a: e) {
//                std::cout << "                         intervals " << a.size() << std::endl;
//            }
//        }
        auto p = get_M_D(v);
        std::cout << "\n E = " << p.first << " D = " << p.second << '\n';
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

BOOST_AUTO_TEST_SUITE_END()
