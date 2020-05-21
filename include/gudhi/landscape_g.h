#include <string>
#include <vector>
#include <limits>
#include <chrono>
#include <random>
#include <atomic>


#include <thread_pool/thread_pool.hpp>
#include <thread>

#include <boost/asio/io_service.hpp>
#include <boost/bind.hpp>
#include <boost/thread/thread.hpp>
#include <boost/asio.hpp>

#include <gudhi/Rips_complex.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>

#include <future>
#include <mutex>

#define DEBUG_GUDHI false

#include "../mean_landscapes/mean_landscapes.h"

namespace smpl {
    using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_fast_persistence>;
    using Filtration_value = Simplex_tree::Filtration_value;
    using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;
    using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
    using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp>;
    using Point = std::vector<double>;


    std::set<int> get_random_sample(int number_of_points,
                               int size_of_one_sample, bool print_pairs = false) {
        std::vector<int> vector_of_points;
        for (int i = 0; i < number_of_points; ++i) {
            vector_of_points.push_back(i);
        }
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(vector_of_points.begin(), vector_of_points.end(), g);

        if (print_pairs) {
            for (const auto &e: vector_of_points) {
                std::cout << e << " ";
            }
            std::cout << std::endl;
        }
        sort(vector_of_points.begin(), vector_of_points.begin() + size_of_one_sample);
//        std::set<int> subcloud(vector_of_points.begin(), vector_of_points.begin() + size_of_one_sample);

        return std::set<int>(vector_of_points.begin(), vector_of_points.begin() + size_of_one_sample);
    }

    std::vector<Point> points_off_reader(const std::string &name_file, double coeff, bool print_pairs = false) {

        std::vector<std::vector<double>> points;
        int number_of_points = 0;
        {
            std::ifstream input_stream(name_file);
            std::string temp_string;
            while (std::getline(input_stream, temp_string)) {
                ++number_of_points;
            }
        }
        if (DEBUG_GUDHI) {
            std::cout << "NUMBER OF POINTS " << number_of_points << std::endl;
        }
        std::ifstream input_stream(name_file);
        std::string line;
        int counter = -1;

    //    std::getline(input_stream, line);

    //    std::cout << line;
    //    std::istringstream s(line);

    //    s >> number_of_points;
        auto subcloud = get_random_sample(number_of_points, (int) (number_of_points * coeff));

        if (print_pairs) {
            std::cout << subcloud.size() << std::endl;
            for (const auto &e: subcloud) {
                std::cout << e << " ";
            }
            std::cout << std::endl;
        }
    //    std::cout << line <<  std::endl;
    //    exit(1);
        while (std::getline(input_stream, line)) {
            ++counter;

            std::vector<double> point;
            std::istringstream s(line);
            if (subcloud.find(counter) == subcloud.end()) {
                continue;
            }
            double value;
            while (s >> value) {
                point.push_back(value);
                s.ignore();
            }
            if (!point.empty()) {
                points.push_back(point);
            }
    //        assert(point.size() == points.front().size());
        }


        if (DEBUG_GUDHI) {
            std::cout << points.size() << std::endl;
            for (const auto &e: points) {

                for (const auto &p: e) {
                    std::cout << p << '_';
                }
                std::cout << points.size() << std::endl;
            }
        }
        if (DEBUG_GUDHI) {
            std::cout << "reducted NUMBER OF POINTS " << points.size() << std::endl;
        }

        return points;
    }


    // infinity
    // Types definition
    std::mutex mute_gudhi;
    void get_diagram(tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> &diagram,
                     std::string filename,
                     int max_rank,
                     double max_edge_length,
                     bool gudhi_format,
                     double subsample_density_coefficient,
                     bool print_pairs) {
        int dim_max = max_rank + 1;
        Rips_complex rips_complex_from_file(points_off_reader(filename, subsample_density_coefficient), max_edge_length,
                                            Gudhi::Euclidean_distance());

        Simplex_tree simplex_tree;
        rips_complex_from_file.create_complex(simplex_tree, dim_max);
        simplex_tree.initialize_filtration();
        mute_gudhi.lock();
        std::cout << simplex_tree.num_simplices() << " simplices" <<  std::endl;
        mute_gudhi.unlock();
        Persistent_cohomology pcoh(simplex_tree);
        pcoh.init_coefficients(2);
        pcoh.compute_persistent_cohomology(1e-10);
        tbb::concurrent_vector<std::vector<std::pair<double, double>>> all_dimension_pairs(dim_max);
        {
    //        pcoh.output_diagram();
            auto &persistence_pairs = pcoh.persistent_pairs_;
            for (const auto &pair: persistence_pairs) {
                int dim = pcoh.cpx_->dimension(get<0>(pair));
                double b = pcoh.cpx_->filtration(get<0>(pair));
                double d = pcoh.cpx_->filtration(get<1>(pair));
                if (d > infinity_1) {
                    d = std::numeric_limits<double>::max();
                }
                all_dimension_pairs[dim].push_back({b, d});
            }

        }
        for (auto &e: all_dimension_pairs) {
            sort(e.begin(), e.end());
        }
        diagram.push_back(all_dimension_pairs);
        if (print_pairs) {
            for (int i = 0; i < all_dimension_pairs.size(); ++i) {
                std::cout << "\nDim " << i <<  std::endl;
                sort(all_dimension_pairs[i].begin(), all_dimension_pairs[i].end());
                for (const auto &e: all_dimension_pairs[i]) {
                    std::cout << e.first << " : " << e.second << std::endl;
                }

            }
        }
    //    exit(17);

    }

//    void program_options(int argc, char *argv[], std::string &off_file_points, std::string &filediag,
//                         Filtration_value &threshold, int &dim_max, int &p, Filtration_value &min_persistence);


    tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>>
    get_diagrams(const std::string &filename,
                 int max_rank,
                 double max_edge_length,
                 bool gudhi_format,
                 int number_of_thread_workers = 1,
                 int number_of_samples = 1,
                 double subsample_density_coefficient = 1.0,
                 bool print_pairs = false) {
    //    tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> diagram;

        tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> all_persistence_diagrams;


    //    boost::asio::thread_pool pool(number_of_thread_workers);
    //    for (int i = 0; i < number_of_samples; ++i) {
    //        boost::asio::post(pool,
    //                          std::bind(get_diagram, std::ref(all_persistence_diagrams), filename, max_rank,
    //                                    max_edge_length, gudhi_format, subsample_density_coefficient, false));
    //    }
    //    pool.join();


    //    tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> all_persistence_diagrams;
        tp::ThreadPoolOptions options;
        options.setThreadCount(number_of_thread_workers);
        tp::ThreadPool pool(options);
        std::vector<std::future<int>> futures(number_of_samples);

        for (int i = 0; i < number_of_samples; ++i) {
            std::packaged_task<int()> t([&filename, &max_rank, &max_edge_length, &gudhi_format, &subsample_density_coefficient, &all_persistence_diagrams]()
                {
    //                get_persistence_pairs_sparse(cloud, radii, subsample_density_coefficient, all_persistence_diagrams);
                    get_diagram(all_persistence_diagrams, filename, max_rank,
                    max_edge_length, gudhi_format, subsample_density_coefficient, false);
                    return 1;
                });
            futures[i] = t.get_future();
            pool.post(t);
        }
        for (int i = 0; i < futures.size(); ++i) {
            int r = futures[i].get();
        }


        std::cout << all_persistence_diagrams.size() << std::endl;
        return all_persistence_diagrams;
    }


    double landscape_gudhi_with_diagrams(std::string from,
                                            std::string to,
                                            tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>>& all_persistence_diagrams,
                                            int max_rank,
                                            double max_edge_length,
                                            int number_of_thread_workers = 1,
                                            int number_of_samples = 1,
                                            double subsample_density_coefficient = 1.0,
                                            bool print_pairs = false,
                                            bool gudhi_format = true) {

        bool print_points = false;
        Filtration_value threshold;
        int dim_max;
        int p;
        Filtration_value min_persistence;
    //    program_options(argc, argv, off_file_points, filediag, threshold, dim_max, p, min_persistence);

        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

        all_persistence_diagrams = get_diagrams(from, max_rank, max_edge_length, true, number_of_thread_workers, number_of_samples, subsample_density_coefficient);
    //    tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> all_persistence_diagrams
    //            = get_diagrams(off_file_points, 2, 0.5, true, 1, 1, 0.4);

        if (print_points) {
            for (int i = 0; i < all_persistence_diagrams.size(); ++i) {
                std::cout << "\n\n\nSAMPLE " << i + 1 <<  std::endl;
                for (int j = 0; j < all_persistence_diagrams[i].size(); ++j) {
                    std::cout << "\nDim " << j <<  std::endl;
    //            sort(all_persistence_diagrams[i][j].begin(), all_persistence_diagrams[i][j].end());
                    for (const auto &e: all_persistence_diagrams[i][j]) {
                        std::cout << e.first << " : " << e.second << std::endl;
                    }

                }
            }
        }

        std::cout << "total samples " << all_persistence_diagrams.size() << std::endl;


        get_average_landscape(all_persistence_diagrams, to);

        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        std::cout << "total duration " << duration <<  std::endl;


        return duration;
    }
double landscape_gudhi(std::string from,
                        std::string to,
                        int max_rank,
                        double max_edge_length,
                        int number_of_thread_workers = 1,
                        int number_of_samples = 1,
                        double subsample_density_coefficient = 1.0,
                        bool print_pairs = false,
                        bool gudhi_format = true) {
        tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> all_persistence_diagrams;
        return landscape_gudhi_with_diagrams(from, to, all_persistence_diagrams, max_rank, max_edge_length,
                number_of_thread_workers, number_of_samples, subsample_density_coefficient, print_pairs, gudhi_format);

    }
}


/*
g++ main.cpp -std=c++17 -lboost_program_options -Os -DNDEBUG -o m
/usr/bin/time -lp  ./m h.off  -r 0.3  -d 3 -p 2


/Users/leonardbee/CLionProjects/Persisitence_Diagram/dots_50


*/