#define BOOST_TEST_MODULE testSuiteCalculator

int num_test_diagram_1 = 10;

#include <atomic>
#include "source/algorithm/landscape_a.cpp"
#include "source/ripser/landscape_r.cpp"
#include <boost/test/included/unit_test.hpp>
#include "../util/compare.cpp"
#include <iostream>

BOOST_AUTO_TEST_SUITE(testSuiteCalculator)

    BOOST_AUTO_TEST_CASE(testOne2) {

            std::string path = "dataset/magnetometer/s2.txt";
            double radii =  0.1;
            //std::cout << "\n\nTest " << path << " diagram with r = 1e11" << std::endl;
            std::vector<double> v;

            tbb::concurrent_vector<std::vector<std::pair<double, double>>> diagram_a;
            tbb::concurrent_vector<std::vector<std::pair<double, double>>> diagram_r;

            for (int i = 0; i < 1; ++i) {
                double time1 = main_ripser(diagram_r, path, "", 2, radii, true, 1, 1, 1, true);
//                //std::cout << " between r  " << diagram_r.size() << std::endl;
//                //std::cout << " before a " << diagram_a.size() << std::endl;
                double time2 = main_algorithm(diagram_a, path, "", 2, radii, true, 1, 1, 1, true);
//                //std::cout << " after r " << diagram_r.size() << std::endl;
//                //std::cout << " after a " << diagram_a.size() << std::endl;

            }

            auto p = get_M_D(v);
//            //std::cout << "r size " << diagram_r.size() << std::endl;
//            for (int i = 0; i < diagram_r.size(); ++i) {
//                //std::cout << " size " << diagram_r[i].size() << std::endl;
//                for (const auto& e: diagram_r[i]) {
//                    //std::cout << e.first << ' ' << e.second << std::endl;
//                }
//            }
//
//            //std::cout << "a size " << diagram_a.size() << std::endl;
//            for (int i = 0; i < diagram_a.size(); ++i) {
//                //std::cout << " size " << diagram_a[i].size() << std::endl;
//                for (const auto& e: diagram_a[i]) {
//                    //std::cout << e.first << ' ' << e.second << std::endl;
//                }
//            }
            auto res = compare(diagram_a, diagram_r);
            BOOST_CHECK_EQUAL(res, true);
            if (!res) {
                exit(1);
            }
    }

    BOOST_AUTO_TEST_CASE(testOne3) {

        std::string path = "dataset/magnetometer/s3.txt";
        double radii = 0.1;
        //std::cout << "\n\nTest " << path << " diagram with r = 1e11" << std::endl;
        std::vector<double> v;

        tbb::concurrent_vector<std::vector<std::pair<double, double>>> diagram_a;
        tbb::concurrent_vector<std::vector<std::pair<double, double>>> diagram_r;

        for (int i = 0; i < 1; ++i) {
//            //std::cout << " before " << diagram_r.size() << std::endl;
//            //std::cout << " before " << diagram_a.size() << std::endl;
            double time1 = main_ripser(diagram_r, path, "", 2, radii, true, 1, 1, 1, true);
//            //std::cout << " between r  " << diagram_r.size() << std::endl;
//            //std::cout << " before a " << diagram_a.size() << std::endl;
            double time2 = main_algorithm(diagram_a, path, "", 2, radii, true, 1, 1, 1, true);
//            //std::cout << " after r " << diagram_r.size() << std::endl;
//            //std::cout << " after a " << diagram_a.size() << std::endl;

        }

        auto p = get_M_D(v);
        //std::cout << "r size " << diagram_r.size() << std::endl;
        for (int i = 0; i < diagram_r.size(); ++i) {
            //std::cout << " size " << diagram_r[i].size() << std::endl;
            for (const auto& e: diagram_r[i]) {
                //std::cout << e.first << ' ' << e.second << std::endl;
            }
        }

        //std::cout << "a size " << diagram_a.size() << std::endl;
        for (int i = 0; i < diagram_a.size(); ++i) {
            //std::cout << " size " << diagram_a[i].size() << std::endl;
            for (const auto& e: diagram_a[i]) {
                //std::cout << e.first << ' ' << e.second << std::endl;
            }
        }
        auto res = compare(diagram_a, diagram_r);
        BOOST_CHECK_EQUAL(res, true);
        if (!res) {
                exit(1);
            }


    }
    

    BOOST_AUTO_TEST_CASE(testEmpty0) {

        std::string path = "dataset/magnetometer/s0.txt";
        double radii = 1e11;
        //std::cout << "\n\nTest " << path << " diagram with r = 1e11" << std::endl;
        std::vector<double> v;

        tbb::concurrent_vector<std::vector<std::pair<double, double>>> diagram_a;
        tbb::concurrent_vector<std::vector<std::pair<double, double>>> diagram_r;

        for (int i = 0; i < 1; ++i) {
//            //std::cout << " before " << diagram_r.size() << std::endl;
//            //std::cout << " before " << diagram_a.size() << std::endl;
            double time1 = main_ripser(diagram_r, path, "", 2, radii, true, 1, 1, 1, true);
//            //std::cout << " betwenn  " << diagram_r.size() << std::endl;
//            //std::cout << " before " << diagram_a.size() << std::endl;
            double time2 = main_algorithm(diagram_a, path, "", 2, radii, true, 1, 1, 1, true);
//            //std::cout << " after  " << diagram_r.size() << std::endl;
//            //std::cout << " before " << diagram_a.size() << std::endl;

        }

        auto p = get_M_D(v);
//        //std::cout << "r size " << diagram_r.size() << std::endl;
//        for (int i = 0; i < diagram_r.size(); ++i) {
//            //std::cout << " size " << diagram_r[i].size() << std::endl;
//            for (const auto& e: diagram_r[i]) {
//                //std::cout << e.first << ' ' << e.second << std::endl;
//            }
//        }
//
//        //std::cout << "a size " << diagram_a.size() << std::endl;
//        for (int i = 0; i < diagram_a.size(); ++i) {
//            //std::cout << " size " << diagram_a[i].size() << std::endl;
//            for (const auto& e: diagram_a[i]) {
//                //std::cout << e.first << ' ' << e.second << std::endl;
//            }
//        }
        auto res = compare(diagram_a, diagram_r);
        BOOST_CHECK_EQUAL(res, true);
        if (!res) {
            exit(1);
        }


    }



   
    BOOST_AUTO_TEST_CASE(testDiagram0) {

            std::string path = "dataset/magnetometer/s50.txt";
            double radii = 1e11;
            //std::cout << "\n\nTest " << path << " diagram with r = 1e11" << std::endl;
            std::vector<double> v;

            tbb::concurrent_vector<std::vector<std::pair<double, double>>> diagram_a;
            tbb::concurrent_vector<std::vector<std::pair<double, double>>> diagram_r;

            for (int i = 0; i < 1; ++i) {
                double time1 = main_ripser(diagram_r, path, "",
                        2, radii, true, 1, 1, 1, true);
                double time2 = main_algorithm(diagram_a, path, "",
                        2, radii, true, 1, 1, 1, true);
            }

            auto p = get_M_D(v);
            auto res = compare(diagram_a, diagram_r);
            BOOST_CHECK_EQUAL(res, true);
            if (!res) {
                exit(1);
            }
    }

    BOOST_AUTO_TEST_CASE(testDiagram1) {

        std::string path = "dataset/figures/dots50.txt";
        double radii = 1e11;
        //std::cout << "\n\nTest " << path << "  diagram with r = 1e11" << std::endl;
        std::vector<double> v;

        tbb::concurrent_vector<std::vector<std::pair<double, double>>> diagram_a;
        tbb::concurrent_vector<std::vector<std::pair<double, double>>> diagram_r;

        for (int i = 0; i < 1; ++i) {

            double time1 = main_ripser(diagram_r, path, "",
                                       2, radii, true, 1, 1, 1, true);
            double time2 = main_algorithm(diagram_a, path, "",
                                          2, radii, true, 1, 1, 1, true);
        }

        auto p = get_M_D(v);
        auto res = compare(diagram_a, diagram_r);

        BOOST_CHECK_EQUAL(res, true);
        if (!res) {
            exit(1);
        }
    }

    BOOST_AUTO_TEST_CASE(testDiagram2) {

        std::string path = "dataset/figures/dots50_no_number.txt";
        double radii = 1e11;
        //std::cout << "\n\nTest " << path << "  diagram with r = 1e11" << std::endl;
        std::vector<double> v;

        tbb::concurrent_vector<std::vector<std::pair<double, double>>> diagram_a;
        tbb::concurrent_vector<std::vector<std::pair<double, double>>> diagram_r;

        for (int i = 0; i < 1; ++i) {

            double time1 = main_ripser(diagram_r, path, "",
                    2, radii, true, 1, 1, 1, true);
            double time2 = main_algorithm(diagram_a, path, "",
                    2, radii, true, 1, 1, 1, true);
        }

        auto p = get_M_D(v);
        auto res = compare(diagram_a, diagram_r);
        BOOST_CHECK_EQUAL(res, true);
        if (!res) {
            exit(1);
        }
    }

    BOOST_AUTO_TEST_CASE(testDiagram3) {

        std::string path = "dataset/magnetometer/s50.txt";
        double radii = 5;
        //std::cout << "\n\nTest " << path << " diagram with r = 5" << std::endl;
        std::vector<double> v;

        tbb::concurrent_vector<std::vector<std::pair<double, double>>> diagram_a;
        tbb::concurrent_vector<std::vector<std::pair<double, double>>> diagram_r;

        for (int i = 0; i < 1; ++i) {

            double time1 = main_ripser(diagram_r, path, "",
                    2, radii, true, 1, 1, 1, true);
            double time2 = main_algorithm(diagram_a, path, "",
                    2, radii, true, 1, 1, 1, true);
        }

        auto p = get_M_D(v);
        auto res = compare(diagram_a, diagram_r);

        BOOST_CHECK_EQUAL(res, true);
        if (!res) {
            exit(1);
        }
    }

    BOOST_AUTO_TEST_CASE(testDiagram4) {

        std::string path = "dataset/figures/dots50.txt";
        double radii = 0.5;
        //std::cout << "\n\nTest " << path << "  diagram with r = 0.5" << std::endl;
        std::vector<double> v;

        tbb::concurrent_vector<std::vector<std::pair<double, double>>> diagram_a;
        tbb::concurrent_vector<std::vector<std::pair<double, double>>> diagram_r;

        for (int i = 0; i < 1; ++i) {

            double time1 = main_ripser(diagram_r, path, "",
                    2, radii, true, 1, 1, 1, true);
            double time2 = main_algorithm(diagram_a, path, "",
                    2, radii, true, 1, 1, 1, true);
        }

        auto p = get_M_D(v);
        auto res = compare(diagram_a, diagram_r);

        BOOST_CHECK_EQUAL(res, true);
        if (!res) {
            exit(1);
        }
    }



   
    BOOST_AUTO_TEST_CASE(testDiagram5) {

            std::string path = "dataset/magnetometer/s1000.txt";
            double radii = 5;
            //std::cout << "\n\nTest " << path << " diagram with r = 5" << std::endl;
            std::vector<double> v;

            tbb::concurrent_vector<std::vector<std::pair<double, double>>> diagram_a;
            tbb::concurrent_vector<std::vector<std::pair<double, double>>> diagram_r;

            for (int i = 0; i < 1; ++i) {

                double time1 = main_ripser(diagram_r, path, "",
                        2, radii, true, 1, 1, 1, true);
                double time2 = main_algorithm(diagram_a, path, "",
                        2, radii, true, 1, 1, 1, true);
            }

            auto p = get_M_D(v);
            auto res = compare(diagram_a, diagram_r);

            BOOST_CHECK_EQUAL(res, true);
            if (!res) {
                exit(1);
            }
    }

    BOOST_AUTO_TEST_CASE(testDiagram6) {

            std::string path = "dataset/figures/bunny500.txt";
            double radii = 0.2;
            //std::cout << "\n\nTest " << path << " diagram with r = 5" << std::endl;
            std::vector<double> v;

            tbb::concurrent_vector<std::vector<std::pair<double, double>>> diagram_a;
            tbb::concurrent_vector<std::vector<std::pair<double, double>>> diagram_r;

            for (int i = 0; i < 1; ++i) {

                double time1 = main_ripser(diagram_r, path, "",
                        2, radii, true, 1, 1, 1, true);
                double time2 = main_algorithm(diagram_a, path, "",
                        2, radii, true, 1, 1, 1, true);
            }

            auto p = get_M_D(v);
            auto res = compare(diagram_a, diagram_r);

            BOOST_CHECK_EQUAL(res, true);
            if (!res) {
                exit(1);
            }
    }

    BOOST_AUTO_TEST_CASE(testDiagram7) {

            std::string path = "dataset/figures/sphere500.txt";
            double radii = 0.2;
            //std::cout << "\n\nTest " << path << " diagram with r = 0.5" << std::endl;
            std::vector<double> v;

            tbb::concurrent_vector<std::vector<std::pair<double, double>>> diagram_a;
            tbb::concurrent_vector<std::vector<std::pair<double, double>>> diagram_r;

            for (int i = 0; i < 1; ++i) {

                double time1 = main_ripser(diagram_r, path, "",
                        2, radii, true, 1, 1, 1, true);
                double time2 = main_algorithm(diagram_a, path, "",
                        2, radii, true, 1, 1, 1, true);
            }

            auto p = get_M_D(v);
            auto res = compare(diagram_a, diagram_r);

            BOOST_CHECK_EQUAL(res, true);
            if (!res) {
                exit(1);
            }
    }

    BOOST_AUTO_TEST_CASE(testDiagram8) {

            std::string path = "dataset/figures/tore500.txt";
            double radii = 0.2;
            //std::cout << "\n\nTest " << path << " diagram with r = 0.5" << std::endl;
            std::vector<double> v;

            tbb::concurrent_vector<std::vector<std::pair<double, double>>> diagram_a;
            tbb::concurrent_vector<std::vector<std::pair<double, double>>> diagram_r;

            for (int i = 0; i < 1; ++i) {

                double time1 = main_ripser(diagram_r, path, "",
                        2, radii, true, 1, 1, 1, true);
                double time2 = main_algorithm(diagram_a, path, "",
                        2, radii, true, 1, 1, 1, true);
            }

            auto p = get_M_D(v);
            auto res = compare(diagram_a, diagram_r);

            BOOST_CHECK_EQUAL(res, true);
            if (!res) {
                exit(1);
            }
    }


    BOOST_AUTO_TEST_CASE(testDiagram9) {

            std::string path = "dataset/figures/human500.txt";
            double radii = 0.2;
            //std::cout << "\n\nTest " << path << " diagram with r = 0.5" << std::endl;
            std::vector<double> v;

            tbb::concurrent_vector<std::vector<std::pair<double, double>>> diagram_a;
            tbb::concurrent_vector<std::vector<std::pair<double, double>>> diagram_r;

            for (int i = 0; i < 5; ++i) {
                double time1 = main_ripser(diagram_r, path, "",
                        2, radii, true, 1, 1, 1, true);
                double time2 = main_algorithm(diagram_a, path, "",
                        2, radii, true, 1, 1, 1, true);
            }

            auto p = get_M_D(v);
//            diagram_a[0][0].first = 188;
            auto res = compare(diagram_a, diagram_r);

            BOOST_CHECK_EQUAL(res, true);
            if (!res) {
                exit(1);
            }
    }


BOOST_AUTO_TEST_SUITE_END()
