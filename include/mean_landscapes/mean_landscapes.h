#pragma once
#include <iostream>
#include <vector>
#include <set>
#include <iterator>
#include <algorithm>
#include <queue>
#include <cmath>
#include <cstring>

#include "tbb/concurrent_vector.h"

#include <gudhi/Persistence_landscape.h>

#define DEBUG_FLAG_MEAN_LANDSCAPES true
#define PLOT_LANDSCAPES true
//using namespace std;

#define infinity_1 (1e30)
#define epsilon_r 1e-10

namespace smpl {
    using Persistence_landscape = Gudhi::Persistence_representations::Persistence_landscape;

    Persistence_landscape compute_average(const std::vector<Persistence_landscape>& to_average) {
        bool dbg = false;

        if (dbg) {
            std::cerr << "to_average.size() : " << to_average.size() << std::endl;
        }

        std::vector<Persistence_landscape> nextLevelMerge(to_average.size());
        for (size_t i = 0; i != to_average.size(); ++i) {
            nextLevelMerge[i] = to_average[i];
        }
        bool is_this_first_level = true;

        while (nextLevelMerge.size() != 1) {
            if (dbg) {
                std::cerr << "nextLevelMerge.size() : " << nextLevelMerge.size() << std::endl;
            }
            std::vector<Persistence_landscape> nextNextLevelMerge;
            nextNextLevelMerge.reserve(to_average.size());
            for (size_t i = 0; i < nextLevelMerge.size(); i = i + 2) {
                if (dbg) {
                    std::cerr << "i : " << i << std::endl;
                }
                Persistence_landscape l;
                if (i + 1 != nextLevelMerge.size()) {
                    l = nextLevelMerge[i] + nextLevelMerge[i + 1];
                } else {
                    l = nextLevelMerge[i];
                }
                nextNextLevelMerge.push_back(l);
            }
            if (dbg) {
                std::cerr << "After this iteration \n";
                getchar();
            }

            is_this_first_level = false;
            nextLevelMerge.swap(nextNextLevelMerge);
        }
        nextLevelMerge[0] *= 1 / static_cast<double>(to_average.size());
        return nextLevelMerge[0];
    }


    std::vector<Persistence_landscape>
    get_average_landscape(
            tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> all_persistence_diagrams,
            std::string path_to_storage = "results",
            bool print_pairs = false) {
        std::cout << "Samples " << all_persistence_diagrams.size() << std::endl;
        for (const auto& e: all_persistence_diagrams) {
            std::cout << "            dims " << e.size() << std::endl;
            for (const auto& a: e) {
                std::cout << "                         intervals " << a.size() << std::endl;
            }
        }
//        exit(1);

        if (DEBUG_FLAG_MEAN_LANDSCAPES) {
            std::cout << "Total size " << all_persistence_diagrams.size() << std::endl;
        }
        if (all_persistence_diagrams.empty()) {
            return std::vector<Persistence_landscape>(0);
        }
//        std::vector<std::vector<Persistence_landscape*>> persistence_landscapes(all_persistence_diagrams[0].size());
        std::vector<std::vector<Persistence_landscape>> persistence_landscapes(all_persistence_diagrams[0].size());

    for (const auto &full_landscape: all_persistence_diagrams) {
            for (int i = 0; i < full_landscape.size(); ++i) {
                print_pairs = true;
                if (DEBUG_FLAG_MEAN_LANDSCAPES) {
                    std::cout << i << " of " << full_landscape.size() << " with " << full_landscape[i].size() << std::endl;
                }
                Persistence_landscape pl;
                if (!full_landscape[i].empty()) {

                    if (i == 0) {
                        pl = Persistence_landscape(full_landscape[i], 3 + 1);
                    } else {
                        //                pl = Persistence_landscape(full_landscape[i], 1); //only upper
                        pl = Persistence_landscape(full_landscape[i], 3);
                    }
                }
                std::string s = "average_" + std::to_string(i);
                const char *ss(s.data());
//                pl.plot(ss);

                persistence_landscapes[i].push_back(pl);
//                persistence_landscapes[i].push_back(&pl);
                if (DEBUG_FLAG_MEAN_LANDSCAPES) {
                    std::cout << "\nDim " << i << std::endl;
                    for (const auto &e: full_landscape[i]) {
                        std::cout << "    " << i << "::: " << e.first << " -> " << e.second << std::endl;
                    }
                }
            }
        }
        std::vector<Persistence_landscape> average_landscape_all_dimensions(all_persistence_diagrams[0].size());
        

        for (int i = 0; i < average_landscape_all_dimensions.size(); ++i) {
            std::vector<Persistence_landscape> plv;
            //std::vector<Persistence_landscape*> plv;
            for (auto &e: persistence_landscapes[i]) {
                plv.push_back(e);
            }
//            average_landscape_all_dimensions[i].compute_average(plv);
            average_landscape_all_dimensions[i] = compute_average(plv);
            std::string s = "average_" + std::to_string(i);
            const char *ss(s.data());

            std::vector<std::string> to(average_landscape_all_dimensions.size());
            for (int i = 0; i < average_landscape_all_dimensions.size(); ++i) {
                std::string new_name = path_to_storage + "/landscape_" + std::to_string(i) + ".txt";
                if (DEBUG_FLAG_MEAN_LANDSCAPES) {
                    std::cout << "to file " << new_name << std::endl;
                }
                {
                    std::ofstream fout(new_name);
                }
                average_landscape_all_dimensions[i].print_to_file(new_name.data());
            }
            if (PLOT_LANDSCAPES) {
                average_landscape_all_dimensions[i].plot(ss);
            }
        }
        return average_landscape_all_dimensions;
    }


    std::pair<double, double> get_M_D(const std::vector<double> &v) {
        double s1 = 0;
        double s2 = 0;

        for (int i = 0; i < v.size(); ++i) {
            s1 += v[i];
            s2 += v[i] * v[i];
        }
        s1 /= v.size();
        s2 = s2 / v.size() - s1 * s1;
        return {s1, sqrt(s2)};
    }
}
