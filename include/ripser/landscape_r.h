#pragma once
#include <iostream>
#include <random>
#include <iterator>
#include <atomic>
#include <boost/asio/io_service.hpp>
#include <boost/bind.hpp>
#include <boost/thread/thread.hpp>
#include <boost/asio.hpp>
#include <algorithm>

//#include <thread_pool/thread_pool.hpp>
#include <thread_pool.hpp>
#include <thread>

#include "ripser.h"
//#include "../mean_landscapes/mean_landscapes.h"

namespace smpl {

    #define DEBUG_RIPSER_0 (false)
//    #define space 5
    //using namespace std;

//    std::atomic<int> r = {0};


    std::set<int> get_random_sample_ripser(int total_size,
                                      int size_of_one_sample, bool print_pairs = false) {
        std::vector<int> vector_of_points(total_size);
        for (int i = 0; i < vector_of_points.size(); ++i) {
            vector_of_points[i] = i;
        }
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(vector_of_points.begin(), vector_of_points.end(), g);
        std::cout << "size of one sample " << size_of_one_sample << std::endl;
        if (DEBUG_RIPSER_0) {
            for (const auto &e: vector_of_points) {
                std::cout << e << " ";
            }
            std::cout << std::endl;
        }
        sort(vector_of_points.begin(), vector_of_points.begin() + size_of_one_sample);
        std::set<int> subcloud(vector_of_points.begin(), vector_of_points.begin() + size_of_one_sample);

        return std::set<int>(vector_of_points.begin(), vector_of_points.begin() + size_of_one_sample);
    }
    tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> get_diagrams_ripser(
            tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>>& diagram,
            const std::string &filename,
            std::string to,
            int max_rank,
            double max_edge_length,
            bool gudhi_format,
            int number_of_thread_workers = 1,
            int number_of_samples = 1,
            double subsample_density_coefficient = 1.0,
            bool print_pairs = true) {
        int number_of_points = 0;

        {
            std::ifstream file_stream(filename);

            std::string line;
            while (std::getline(file_stream, line)) {
                ++number_of_points;
                std::vector<value_t> point;
                std::istringstream s(line);
            }
        }
        if (number_of_points == 0) {
            tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> tbcv (1, tbb::concurrent_vector<std::vector<std::pair<double, double>>>(max_rank + 1));
            diagram = tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> (max_rank + 1);
            return tbcv;
        }
    //    if (gudhi_format && number_of_points > 0) {
    //        --number_of_points;
    //    }
        if (DEBUG_RIPSER_0) {
            std::cout << "total number of points " << number_of_points << std::endl;
        }

        std::vector<int> cloud(number_of_points);
        for (int i = 0; i < number_of_points; ++i) {
            cloud[i] = i;
        }
        std::vector<std::string> argv_strings = {"./ripser", "--dim", std::to_string(max_rank), "--threshold",
                                                  std::to_string(max_edge_length),
                                                  "--modulus", "2", "--format", "point-cloud",
                                                  filename};
    //    int number_of_thread_workers = 1;
    //
    //    double subsample_density_coefficient = 1;
    //    int number_of_samples = 1;
    //    boost::asio::thread_pool pool(number_of_thread_workers);
    //    tbb::concurrent_vector < tbb::concurrent_vector < std::vector < std::pair < double, double >>
    //                                                                                               >> all_persistence_diagrams;
    //    for (int i = 0; i < number_of_samples; ++i) {
    //        auto sample = get_random_sample_ripser(cloud, (int) (cloud.size() * subsample_density_coefficient));
    //        boost::asio::post(pool,
    //                          std::bind(main_ripser_init, 10, argv_strings, sample,
    //                                    std::ref(all_persistence_diagrams)));
    //    }
    //    pool.join();
        tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> all_persistence_diagrams;
        tp::ThreadPoolOptions options;
        options.setThreadCount(number_of_thread_workers);
        tp::ThreadPool pool(options);
        std::vector<std::future<int>> futures(number_of_samples);

        for (int i = 0; i < number_of_samples; ++i) {
            std::packaged_task<int()> t([&number_of_points, &argv_strings, &all_persistence_diagrams, &subsample_density_coefficient]()
                {
                    auto sample = get_random_sample_ripser(number_of_points, (int) (number_of_points * subsample_density_coefficient));
                    std::cout << "sample size " <<  sample.size() << std::endl;
                    main_ripser_init(10, argv_strings, sample, all_persistence_diagrams);
                    return 1;
                });

            futures[i] = t.get_future();
            pool.post(t);
        }
        for (int i = 0; i < futures.size(); ++i) {
            int r = futures[i].get();
        }




        if (DEBUG_RIPSER_0) {
            std::cout << number_of_thread_workers << " on samples: " << number_of_samples << std::endl;
            std::cout << "____________________________________________________\n\n\n\n" << all_persistence_diagrams.size()
                 << std::endl;
        }

        for (int i = 0; i < all_persistence_diagrams.size(); ++i) {
            if (DEBUG_RIPSER_0) {
                std::cout << "\n\n\n\ncurrently in " << i << std::endl;
            }
            for (int j = 0; j < all_persistence_diagrams[i].size(); ++j) {
                if (DEBUG_RIPSER_0) {
                    std::cout << "                 into j " << j << std::endl;
                }
                for (auto &e: all_persistence_diagrams[i][j]) {

    //                if (e.second > 1e10) {
    //                    e.second = std::numeric_limits<double>::max();
    //                }
                    sort(all_persistence_diagrams[i][j].begin(), all_persistence_diagrams[i][j].end());
                    if (DEBUG_RIPSER_0) {
                        std::cout << e.first << " and " << e.second << std::endl;
                    }
                }
            }

        }
        diagram = all_persistence_diagrams;
//        get_average_landscape(all_persistence_diagrams, to);
//        if (!all_persistence_diagrams.empty()) {
//            diagram = all_persistence_diagrams;
//        }
        return all_persistence_diagrams;
    }


    double main_ripser(tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>>& diagram,
                        std::string from,
                        std::string to,
                        int max_rank,
                        double max_edge_length,
                        int number_of_thread_workers = 1,
                        int number_of_samples = 1,
                        double subsample_density_coefficient = 1.0,
                        bool print_pairs = false,
                        bool gudhi_format = true) {
        std::vector<int> v;
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        if (subsample_density_coefficient < 0 || subsample_density_coefficient > 1) {
            subsample_density_coefficient = 1;
        }
    //    get_diagrams_ripser(filetore, 2, 0.5, true, 10, 10, 0.4, true); //not sampled
        get_diagrams_ripser(diagram, from, to, max_rank, max_edge_length, true, number_of_thread_workers, number_of_samples,
                     subsample_density_coefficient, print_pairs);
    //    while (diagram.size() < max_rank + 1) {
    //        diagram.push_back(std::vector<std::pair<double, double>> (0));
    //    }
        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        std::cout << "total duration " << duration << std::endl;
        return duration;
    }
    double landscape_ripser(std::string from,
                        std::string to,
                        int max_rank,
                        double max_edge_length,
                        int number_of_thread_workers = 1,
                        int number_of_samples = 1,
                        double subsample_density_coefficient = 1.0) {
        tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> tmp;
        return main_ripser(tmp, from, to, max_rank, max_edge_length, number_of_thread_workers, number_of_samples, subsample_density_coefficient, true, true);
    }
}