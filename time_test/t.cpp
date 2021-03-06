#define BOOST_TEST_MODULE testSuiteCalculator

int num_ = 10;
#include <boost/test/included/unit_test.hpp>
//#include "/Users/leonardbee/CLionProjects/SubsamplingMethodsForPersistenceLandscape1/tests.cpp"
#include <iostream>

using namespace smpl;

BOOST_AUTO_TEST_SUITE(testSuiteCalculator)
    
    BOOST_AUTO_TEST_CASE(testCalculator1) {
        std::ofstream out;

        out.open("log.txt", std::ios::app);
        out << "\n\n\nnew log\n r=0.5 sampled ";
        std::cout << "\n\nTest sampled diagram with r = 0.5" << std::endl;
        std::vector<double> v;
        for (int i = 0; i < num_; ++i) {
            double time = main_ripser("..dataset/human500.txt", "/Users/leonardbee/Desktop/dataset/tore/sampled_persistence",
                        2, 0.5, true, 4, 10, 0.4);
            v.push_back(time);
            out << time << ' ';
        }
        auto p = get_M_D(v);
        out << "\n E = " << p.first << " D = " << p.second << '\n';
        std::cout << "\n E = " << p.first << " D = " << p.second << '\n';
        out.close();
        
    }
    BOOST_AUTO_TEST_CASE(testCalculator2) {
        std::ofstream out;
        out.open("log.txt", std::ios::app);
        out << "\n\n\nnew log\n r=0.5 full ";
        std::cout << "\n\nTest full diagram with r = 0.5" << std::endl;
        std::vector<double> v;
        for (int i = 0; i < num_; ++i) {
            double time = main_ripser("..dataset/human500.txt", "/Users/leonardbee/Desktop/dataset/tore/sampled_persistence",
                        2, 0.5, true, 1, 1, 1);
            v.push_back(time);
            out << time << ' ';
        }
        auto p = get_M_D(v);
        out << "\n E = " << p.first << " D = " << p.second << '\n';
        std::cout << "\n E = " << p.first << " D = " << p.second << '\n';
        out.close();
    }

    
    BOOST_AUTO_TEST_CASE(testCalculator3) {
        std::ofstream out;
        out.open("log.txt", std::ios::app);
        out << "\n\n\nnew log\n r=diam sampled ";
        std::cout << "\n\nTest sampled diagram with r =  diam" << std::endl;
        std::vector<double> v;
        for (int i = 0; i < num_; ++i) {
            double time = main_ripser("..dataset/human500.txt", "/Users/leonardbee/Desktop/dataset/tore/sampled_persistence",
                        2, 1e11, true, 4, 10, 0.4);
            v.push_back(time);
            out << time << ' ';
        }
        auto p = get_M_D(v);
        out << "\n E = " << p.first << " D = " << p.second << '\n';
        std::cout << "\n E = " << p.first << " D = " << p.second << '\n';
        out.close();
    }
    BOOST_AUTO_TEST_CASE(testCalculator4) {
        std::ofstream out;
        out.open("log.txt", std::ios::app);
        out << "\n\n\nnew log\n r=diam full ";
        std::cout << "\n\nTest full diagram with r = 0.5" << std::endl;
        std::vector<double> v;
        for (int i = 0; i < num_; ++i) {
            double time = main_ripser("..dataset/human500.txt", "/Users/leonardbee/Desktop/dataset/tore/sampled_persistence",
                        2, 1e11, true, 1, 1, 1);
            v.push_back(time);
            out << time << ' ';
        }
        auto p = get_M_D(v);
        out << "\n E = " << p.first << " D = " << p.second << '\n';
        std::cout << "\n E = " << p.first << " D = " << p.second << '\n';
        out.close();
//        BOOST_CHECK_EQUAL(s(6,6), 12);
    }

BOOST_AUTO_TEST_SUITE_END()
