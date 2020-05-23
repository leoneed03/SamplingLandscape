#pragma once
#ifndef TRAVIS_OK_LANDSCAPES_ALGORITHM_H
#define TRAVIS_OK_LANDSCAPES_ALGORITHM_H

#include "../mean_landscapes/mean_landscapes.h"
#include "landscape_a.h"

namespace smpl {
    double landscape_with_algorithm(std::string from,
                                    std::string to,
                                    int max_rank,
                                    double max_edge_length,
                                    int number_of_thread_workers = 1,
                                    int number_of_samples = 1,
                                    double subsample_density_coefficient = 1.0,
                                    bool print_pairs = false,
                                    bool gudhi_format = true) {

        tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> apd;
        double time = landscape_algorithm_with_diagrams(from, to, apd, max_rank, max_edge_length,
                                                        number_of_thread_workers,
                                                        number_of_samples, subsample_density_coefficient);
        get_average_landscape(apd, to);
        return time;
    }
}
#endif //TRAVIS_OK_LANDSCAPES_ALGORITHM_H
