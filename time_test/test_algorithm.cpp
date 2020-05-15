#define BOOST_TEST_MODULE testSuiteCalculator

int num_ = 10;

#include <boost/test/included/unit_test.hpp>
#include "source/algorithm/landscape_a.cpp"
#include <iostream>
BOOST_AUTO_TEST_SUITE(testSuiteCalculator)

    BOOST_AUTO_TEST_CASE(testCalculator0) {
            std::ofstream out;
            out.open("log.txt", std::ios::app);
            out << "\n\n\nnew log\n r=0.5 sampled ";
            std::cout << "\n\nTest sampled diagram with r = 0.5" << std::endl;
            std::vector<double> v;
            for (int i = 0; i < 1; ++i) {
    //            double time = main_ripser("dataset/figures/dots50_no_number.txt", "/Users/leonardbee/Desktop/dataset/tore/sampled_persistence",
    //                        2, 8, true, 1, 1, 1, true);
                double time = main_algorithm("dataset/magnetometer/s50.txt", "/Users/leonardbee/Desktop/dataset/tore/sampled_persistence",
                                          2, 1e11, true, 1, 1, 1, true);
                v.push_back(time);
                out << time << ' ';
            }
            auto p = get_M_D(v);
            out << "\n E = " << p.first << " D = " << p.second << '\n';
            std::cout << "\n E = " << p.first << " D = " << p.second << '\n';
            out.close();

    }



    BOOST_AUTO_TEST_CASE(testCalculator1) {
            std::ofstream out;
            std::cout << "\n\nTest sampled diagram with r = 8" << std::endl;
            std::vector<double> v, av;
            for (int i = 0; i <  1 ; ++i) {
                std::clog << "\n\nSTARTED " << i << std::endl;
                double time = main_algorithm("dataset/figures/dots50_no_number.txt", "",
                                             2, 8, true, 1,1,1);

                std::clog << "\n\nCOMPUTED " << i << std::endl;
                v.push_back(time);
                av.push_back(matrix_size_cntr);
                out << time << ' ';
            }
            auto p = get_M_D(v);
            auto total_size = get_M_D(av);
            out << "\n E = " << p.first << " D = " << p.second << '\n';
            std::cout << "\n E = " << p.first << " D = " << p.second << '\n';
            std::cout << "\n size E = " << total_size.first << " size D = " << total_size.second << '\n';
            out.close();

    }


    /*
    BOOST_AUTO_TEST_CASE(testCalculator1) {
        std::ofstream out;
        out.open("logalgorithm.txt", std::ios::app);
        out << "\n\n\nnew logalgorithm.n r=0.5 sampled ";
        std::cout << "\n\nTest sampled diagram with r = 0.5" << std::endl;
        std::vector<double> v, av;
        for (int i = 0; i < num_ ; ++i) {
               std::clog << "\n\nSTARTED " << i << std::endl;
            double time = main_algorithm("dataset/figures/human500.txt", "",
                        2, 0.5, true, 4, 7, 0.1);
            
            std::clog << "\n\nCOMPUTED " << i << std::endl;
            v.push_back(time);
            av.push_back(matrix_size_cntr);
            out << time << ' ';
        }
        auto p = get_M_D(v);
        auto total_size = get_M_D(av);
        out << "\n E = " << p.first << " D = " << p.second << '\n';
        std::cout << "\n E = " << p.first << " D = " << p.second << '\n';
        std::cout << "\n size E = " << total_size.first << " size D = " << total_size.second << '\n';
        out.close();
        
    }*/

BOOST_AUTO_TEST_SUITE_END()
