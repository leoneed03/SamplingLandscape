#define BOOST_TEST_MODULE testSuiteCalculator

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <random>
#include <iterator>
#include <atomic>
#include <boost/asio/io_service.hpp>
#include <boost/bind.hpp>
#include <boost/thread/thread.hpp>
#include <boost/asio.hpp>
#include <algorithm>

#include "ripser.cpp"
#include "../mean_landscapes/mean_landscapes.cpp"

#define space 5
using namespace std;

std::atomic<int> r = 0;


set<int> get_random_sample_ripser(vector<int> &vector_of_points,
                                  int size_of_one_sample, bool print_pairs = false) {
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(vector_of_points.begin(), vector_of_points.end(), g);

    if (print_pairs) {
        for (const auto &e: vector_of_points) {
            cout << e << " ";
        }
        cout << endl;
    }
    sort(vector_of_points.begin(), vector_of_points.begin() + size_of_one_sample);
    std::set<int> subcloud(vector_of_points.begin(), vector_of_points.begin() + size_of_one_sample);

    return std::set<int>(vector_of_points.begin(), vector_of_points.begin() + size_of_one_sample);
}

tbb::concurrent_vector <tbb::concurrent_vector<std::vector < std::pair < double, double>>>>

get_diagrams_ripser(tbb::concurrent_vector<std::vector<std::pair<double, double>>>& diagram,
                    const std::string &filename,
                    int max_rank,
                    double max_edge_length,
                    bool gudhi_format,
                    int number_of_thread_workers = 1,
                    int number_of_samples = 1,
                    double subsample_density_coefficient = 1.0,
                    bool print_pairs = false) {

    std::ifstream file_stream(filename);
    int number_of_points = 0;
    std::string line;
    while (std::getline(file_stream, line)) {
        ++number_of_points;
        std::vector <value_t> point;
        std::istringstream s(line);
    }
//    if (gudhi_format && number_of_points > 0) {
//        --number_of_points;
//    }
    if (print_pairs) {
        std::cout << "total number of points " << number_of_points << std::endl;
    }

    vector<int> cloud(number_of_points);
    for (int i = 0; i < number_of_points; ++i) {
        cloud[i] = i;
    }
    std::vector <std::string> argv_strings = {"./ripser", "--dim", std::to_string(max_rank), "--threshold",
                                              std::to_string(max_edge_length),
                                              "--modulus", "2", "--format", "point-cloud",
                                              filename};
//    int number_of_thread_workers = 1;
//
//    double subsample_density_coefficient = 1;
//    int number_of_samples = 1;
    boost::asio::thread_pool pool(number_of_thread_workers);
    tbb::concurrent_vector < tbb::concurrent_vector < std::vector < std::pair < double, double >>
                                                                                               >> all_persistence_diagrams;
    for (int i = 0; i < number_of_samples; ++i) {
        auto sample = get_random_sample_ripser(cloud, (int) (cloud.size() * subsample_density_coefficient));
        boost::asio::post(pool,
                          std::bind(main_ripser_init, 10, argv_strings, sample,
                                    std::ref(all_persistence_diagrams)));
    }
    pool.join();
    cout << number_of_thread_workers << " on samples: " << number_of_samples << endl;
    if (print_pairs) {
        cout << "____________________________________________________\n\n\n\n" << all_persistence_diagrams.size()
             << endl;
    }
    for (int i = 0; i < all_persistence_diagrams.size(); ++i) {
        diagram = all_persistence_diagrams[i];
        if (print_pairs) {
            cout << "\n\n\n\ncurrently in " << i << endl;
        }
        for (int j = 0; j < all_persistence_diagrams[i].size(); ++j) {
            if (print_pairs) {
                cout << "                 into j " << j << endl;
            }
            for (auto &e: all_persistence_diagrams[i][j]) {

//                if (e.second > 1e10) {
//                    e.second = std::numeric_limits<double>::max();
//                }
                sort(all_persistence_diagrams[i][j].begin(), all_persistence_diagrams[i][j].end());
                if (print_pairs) {
                    cout << e.first << " and " << e.second << endl;
                }
            }
        }

    }
    get_average_landscape(all_persistence_diagrams, "/Users/leonardbee/Desktop/dataset/tore/sampled_persistence");
    if (!all_persistence_diagrams.empty()) {
        diagram = all_persistence_diagrams[0];
    }
    return all_persistence_diagrams;
}

double main_ripser(tbb::concurrent_vector<std::vector<std::pair<double, double>>>& diagram,
                    std::string from, std::string to,
                    int max_rank,
                    double max_edge_length,
                    bool gudhi_format,
                    int number_of_thread_workers = 1,
                    int number_of_samples = 1,
                    double subsample_density_coefficient = 1.0,
                    bool print_pairs = false) {
    vector<int> v;
//
//    std::string filetore = "/Users/leonardbee/Desktop/dataset/human500.txt";
//
//    std::string filetore = "/Users/leonardbee/CLionProjects/SubsamplingMethodsForPersistenceLandscape1/cloud_50";
//    std::string filename = "/Users/leonardbee/CLionProjects/SubsamplingMethodsForPersistenceLandscape1/dots_50.txt";

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
//    get_diagrams_ripser(filetore, 2, 0.5, true, 10, 10, 0.4, true); //not sampled
    get_diagrams_ripser(diagram, from, max_rank, max_edge_length, true, number_of_thread_workers, number_of_samples,
                 subsample_density_coefficient, print_pairs);
//    get_diagrams_ripser(filetore, 2, 0.5, true); //sampled  persistence
//    get_diagrams_ripser(filetore, 2, 0.5, true, 4, 10, 0.4); //sampled  persistence
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << "total duration " << duration << std::endl;
    return duration;
}