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
#include "mean_landscapes.cpp"

#define space 5
using namespace std;

std::atomic<int> r = 0;

void p() {
    std::cout << "Hello there " << std::endl;
}

set<int> get_random_sample(vector<int> &vector_of_points, int size_of_one_sample) {
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(vector_of_points.begin(), vector_of_points.end(), g);
    for (const auto &e: vector_of_points) {
        cout << e << " ";
    }
    cout << endl;
    sort(vector_of_points.begin(), vector_of_points.begin() + size_of_one_sample);
//    ripser_func();
    std::set<int> subcloud(vector_of_points.begin(), vector_of_points.begin() + size_of_one_sample);

    return std::set<int>(vector_of_points.begin(), vector_of_points.begin() + size_of_one_sample);
}


int main() {
    vector<int> v;
    int n = 10;

    for (int i = 0; i < n; ++i) {
//        v.emplace_back(i);
    }
    for (int i = 0; i < n; ++i) {
//        get_random_sample(v);
    }


    double max_edge_length = 8;
    double min_edge_length = 1e-10;
    int max_cohomology_rank = 2;


    int max_rank = 3;


    std::string filename = "/Users/leonardbee/CLionProjects/restribution/cloud_50";
    std::ifstream file_stream(filename);
    int number_of_points = 0;
    std::string line;
    while (std::getline(file_stream, line)) {
        ++number_of_points;
        std::vector<value_t> point;
        std::istringstream s(line);
    }
    std::cout << "total number of points " << number_of_points << std::endl;

    vector<int> cloud(number_of_points);
    for (int i = 0; i < number_of_points; ++i) {
        cloud[i] = i;
    }
    std::vector<std::string> argv_strings = {"./ripser", "--dim", std::to_string(max_rank), "--threshold", std::to_string(max_edge_length),
                                             "--modulus", "2", "--format", "point-cloud",
                                             "/Users/leonardbee/CLionProjects/restribution/cloud_50"};
//    auto res = main_ripser(10, argv_strings, get_random_sample(cloud, (int) (cloud.size() * 1)));
//    std::cout << res.size();
//    for (const auto &dim: res) {
//        for (const auto &persistence_pair: dim) {
//            std::cout << persistence_pair.first << " vs " << persistence_pair.second << std::endl;
//        }
//        std::cout << "_______\n_______\n____________________________________\n" << std::endl;
//    }






//
//    boost::asio::io_service ioService;
//    boost::thread_group threadpool;
//
//    boost::asio::io_service::work work(ioService);
    int number_of_thread_workers = 4;

    double subsample_density_coefficient = 0.3;
    int number_of_samples = 30;
    boost::asio::thread_pool pool(number_of_thread_workers);
    tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> all_persistence_diagrams(
            number_of_thread_workers);
    for (int i = 0; i < number_of_samples; ++i) {
        auto sample = get_random_sample(cloud, (int) (cloud.size() * subsample_density_coefficient));
        boost::asio::post(pool,
                          std::bind(main_ripser, 10, argv_strings, sample,
                                    std::ref(all_persistence_diagrams)));
//        for (int i = 0; i < current_persistence_diagrams.size(); ++i) {
//            all_persistence_diagrams[i].emplace_back(current_persistence_diagrams[i]);
//        }
    }



//    auto resu = get_random_sample(cloud, (int) (cloud.size() * 1));
//    boost::asio::post(pool, std::bind(main_ripser, 10, argv_strings, std::ref(resu)));
//    boost::asio::post(pool, std::bind(main_ripser, 10, argv_strings, std::ref(resu)));
//    boost::asio::post(pool, p);
    pool.join();
    std::cout << "counter " << special_counter << std::endl;
    for (int i = 0; i < all_persistence_diagrams.size(); ++i) {
        std::cout << "dimension " << i << " ";

        if (all_persistence_diagrams[i].size() != number_of_samples) {
            std::cout << all_persistence_diagrams[i].size() << std::endl;
            exit(19);
        }
        for (auto &vector_of_intervals: all_persistence_diagrams[i]) {
            for (auto& persistence_interval: vector_of_intervals) {
                std::cout << endl << persistence_interval.first << " -> before -> " << persistence_interval.second << endl;
                std::pair<double, double> new_pair = {(persistence_interval.second - persistence_interval.first) / 2.0, (persistence_interval.second + persistence_interval.first) / 2.0};


                std::swap(new_pair, persistence_interval);
                std::cout << endl << persistence_interval.first << " -> after -> " << persistence_interval.second << endl;
            }

        }
        cout << endl;
    }

    std::vector<std::vector<std::pair<double, double>>> mean_landscape;
    for (auto& landscape: all_persistence_diagrams) {
        mean_landscape.emplace_back(get_mean_persistence(landscape));
    }

    for (int i = 0; i < mean_landscape.size(); ++i) {
        cout << i << ": " << mean_landscape[i].size() << endl;
    }




    return 0;
}