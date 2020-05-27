#pragma once

#include <phat/compute_persistence_pairs.h>
#include <phat/algorithms/chunk_reduction.h>

#include <future>
#include <atomic>
#include <iostream>
#include <vector>
#include <algorithm>
#include <thread>
#include <set>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <map>
#include <thread>
#include <shared_mutex>
#include <mutex>
#include <queue>

#include <thread_pool.hpp>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/concurrent_vector.h>
#include <boost/container/flat_map.hpp>
#include <unordered_map>
#include <chrono>
#include <random>


#include <boost/bind.hpp>
#include <boost/thread/thread.hpp>
#include <boost/asio.hpp>
#include <boost/dynamic_bitset.hpp>

#include "Subcloud.h"

namespace smpl {

    std::ostream &operator<<(std::ostream &os, const std::vector<std::vector<double>> &matrix) {
        for (const auto &row: matrix) {
            for (const auto &element: row) {
                os << element << " ";
            }
            os << std::endl;
        }
        return os;
    }


    std::vector<int> get_random_sample(std::vector<int> &vector_of_points, int size_of_one_sample) {
        std::random_device rd;
        std::mt19937 g(rd());
        shuffle(vector_of_points.begin(), vector_of_points.end(), g);
        sort(vector_of_points.begin(), vector_of_points.begin() + size_of_one_sample);
        std::vector<int> subcloud(vector_of_points.begin(), vector_of_points.begin() + size_of_one_sample);
        std::vector<int> set_of_points = std::vector<int>(vector_of_points.begin(),
                                                          vector_of_points.begin() + size_of_one_sample);
        return set_of_points;
    }

    void get_persistence_pairs_sparse(Cloud *cloud, double radii, double subsample_density_coefficient,
                                      tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> &v_pairs) {

        if (DEBUG_FLAG_0) {
            //mute.lock();
            std::cout << "non flag entered pairs construction" << std::endl;
            std::cout << "entered pairs construction" << std::endl;
            //mute.unlock();
        }
        int max_number_of_points_in_simplex = 4;
        int number_of_dots = cloud->get_size();
        std::vector<int> dots(number_of_dots);

        for (int i = 0; i < number_of_dots; ++i) {
            dots[i] = i;
        }

        auto subsample = get_random_sample(dots, (int) ((dots.size()) * subsample_density_coefficient));
        if (DEBUG_FLAG_0) {
            //mute.lock();
            std::cout << "trying construct tree" << std::endl;
            //mute.unlock();
        }
        auto subcloud = new SubCloud<phat::chunk_reduction, phat::vector_vector>(*cloud, subsample, radii, max_number_of_points_in_simplex);

        if (DEBUG_FLAG_0) {
            //mute.lock();
            std::cout << "constructed tree" << std::endl;
            //mute.unlock();
        }


        if (DEBUG_FLAG_0) {
            std::cout << number_of_dots << " started building tree \n\n\n\n" << subsample.size() << std::endl;
            for (const auto &e: subsample) {
                std::cout << e << ',';
            }
            std::cout << "_________________________________________________________________________" << std::endl;
        }


        subcloud->insert_all_simplices_including(max_number_of_points_in_simplex);
        std::vector<std::vector<std::pair<double, double>>> result;

        //mute.lock();

        int ss = 0;
        for (const auto &e: (subcloud->get_root()->get_simplices()  )) {
            ss += e.size();
        }

        if (DEBUG_FLAG_0) {
            std::cout << "started calculating boundary matrix compressed " << ss << std::endl;
        }
//        matrix_size_cntr += ss;
        //mute.unlock();
        //only  1


        auto Persistence_landscape_a = subcloud->get_all_dimensions_landscape(); //DEBUG_FLAG_0 cohomology //58 Mb memory 332 (282) ms (50 dots r = 8)

        if (DEBUG_FLAG_0) {
            std::cout << "delete subcloud root" << std::endl;
            std::cout << "delete subcloud" << std::endl;
        }
        delete subcloud;

        //mute.lock();
        for (int i = 0; i < Persistence_landscape_a.size(); ++i) {
            if (DEBUG_FLAG_0) {
                std::cout << "DIM_" << i << std::endl;
            }
            sort(Persistence_landscape_a[i].begin(), Persistence_landscape_a[i].end());
            if (DEBUG_FLAG_0) {
                for (const auto &e: Persistence_landscape_a[i]) {
                    std::cout << e.first << " " << e.second << std::endl;
                }
            }
        }
        v_pairs.push_back(Persistence_landscape_a);
        //mute.unlock();
    }

    void get_average_landscape_once(std::string path,
                                    tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> &diagrams,
                                    Cloud *cloud,
                                    int number_of_thread_workers,
                                    double radii = 0.5,
                                    double subsample_density_coefficient = 0.3,
                                    int number_of_samples = 10,
                                    bool print_pairs = false) {

        //    boost::asio::thread_pool pool(number_of_thread_workers);
        if (DEBUG_FLAG_0) {
            std::cout << "Created pool" << std::endl;
        }
        //    for (int i = 0; i < number_of_samples; ++i) {
        //        boost::asio::post(pool,
        //                          bind(get_persistence_pairs_sparse, cloud, radii, subsample_density_coefficient, ref(all_persistence_diagrams)));
        //    }
        //    pool.join();




        tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> all_persistence_diagrams;
        tp::ThreadPoolOptions options;
        options.setThreadCount(number_of_thread_workers);
        tp::ThreadPool pool(options);
        std::vector<std::future<int>> futures(number_of_samples);

        for (int i = 0; i < number_of_samples; ++i) {
            std::packaged_task<int()> t([&cloud, &radii, &subsample_density_coefficient, &all_persistence_diagrams]() {
                get_persistence_pairs_sparse(cloud, radii, subsample_density_coefficient, all_persistence_diagrams);
                return 1;
            });
            futures[i] = t.get_future();
            pool.post(t);
        }
        for (int i = 0; i < futures.size(); ++i) {
            int r = futures[i].get();
        }

        {
            if (DEBUG_FLAG_0) {
                //mute.lock();
                std::cout << "joined!" << std::endl;
                //mute.unlock();
            }
        }
        int counter = 1;
        if (all_persistence_diagrams.empty()) {
            if (DEBUG_FLAG_0) {
                std::cout << "no landscapes" << std::endl;
            }
            return;
        }
        diagrams = all_persistence_diagrams;
//        get_average_landscape(all_persistence_diagrams, path);
    }

    double
    main_algorithm(tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> &diagrams,
                   std::string from,
                   std::string to,
                   int max_rank,
                   double max_edge_length,
                   int number_of_thread_workers = 1,
                   int number_of_samples = 1,
                   double subsample_density_coefficient = 1.0,
                   bool print_pairs = false,
                   bool gudhi_format = true) {
//        zero_cntr = 0;
//        extra_cntr = 0;
//        matrix_size_cntr = 0;
        int max_number_of_points_in_simplex = max_rank + 2;

        if (DEBUG_FLAG_0) {
            std::cout << "Trying to create" << std::endl;
        }

        Cloud *matrix = new Cloud(from, max_number_of_points_in_simplex);
        if (DEBUG_FLAG_0) {
            std::cout << "Created" << std::endl;
        }
        std::vector<int> dots;

        int number_of_dots = matrix->get_size();
        std::set<int> not_include = {0, 2, 6, 8, 23, 25};
        for (int i = 0; i < number_of_dots; ++i) {
            if (not_include.find(i) == not_include.end()) {
                dots.emplace_back(i);
            }
        }

        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        if (DEBUG_FLAG_0) {
            std::cout << "Started calculating" << std::endl;
        }
        get_average_landscape_once(to, diagrams, matrix, number_of_thread_workers, max_edge_length,
                                   subsample_density_coefficient,
                                   number_of_samples);
        if (DEBUG_FLAG_0) {
            std::cout << "Finished calculating" << std::endl;
        }
        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();


        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

        if (DEBUG_FLAG_0) {
            std::cout << "Deleting tree" << std::endl;
        }
        delete matrix;

        //mute.lock();
//        std::cout << zero_cntr << " so " << extra_cntr << std::endl << "duration " << duration << std::endl;
        //mute.unlock();

        return duration;
    }

    double landscape_algorithm(std::string from,
                               std::string to,
                               int max_rank,
                               double max_edge_length,
                               int number_of_thread_workers,
                               int number_of_samples,
                               double subsample_density_coefficient,
                               bool print_pairs = false,
                               bool gudhi_format = true) {
        tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> tmp;
        return main_algorithm(tmp, from, to, max_rank, max_edge_length,
                              number_of_thread_workers, number_of_samples, subsample_density_coefficient, false, true);

    }

    double landscape_algorithm_with_diagrams(std::string from,
                                             std::string to,
                                             tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> &all_persistence_diagrams,
                                             int max_rank,
                                             double max_edge_length,
                                             int number_of_thread_workers,
                                             int number_of_samples,
                                             double subsample_density_coefficient,
                                             bool print_pairs = false,
                                             bool gudhi_format = true) {
        return main_algorithm(all_persistence_diagrams, from, to, max_rank, max_edge_length,
                              number_of_thread_workers, number_of_samples, subsample_density_coefficient, false, true);
    }
}

