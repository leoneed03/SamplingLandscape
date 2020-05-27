#pragma once
#ifndef CLOUD_H
#define CLOUD_H

#include "Simplex_tree_owner.h"

namespace smpl {
    struct Cloud {
    private:
        int dimension, size;
        std::vector<std::vector<double>> points;
        std::vector<std::vector<double>> distances;
        Simplex_tree_owner *simplex_tree;

    public:
        int get_size() {
            return size;
        }
        const std::vector<std::vector<double>>& get_distances() const {
            return distances;
        }
        Simplex_tree_owner * get_simplex_tree_owner() const {
            return simplex_tree;
        }
        Cloud(const std::string &path, int n) {
            std::ifstream input_stream(path);
            std::cout << path << std::endl;
            if (!input_stream) {
                std::cout << "problems opening file" << std::endl;
                exit(1);
            }

            int counter = 0;
            int dim = -1;
            std::string line;
            std::vector<std::vector<double>> points;

            while (std::getline(input_stream, line)) {
                ++counter;
                std::vector<double> point;
                std::istringstream s(line);
                double value;
                while (s >> value) {
                    point.push_back(value);
                    s.ignore();
                }

                if (!point.empty()) {
                    if (dim == -1) {
                        dim = point.size();
                        points.push_back(point);
                        continue;
                    }
                    if (dim != point.size()) {
                        std::cout << "Wrong file format" << std::endl;
                        exit(10);
                    }
                    points.push_back(point);
                }
                //            assert(point.size() == points.front().size());
            }

            size = counter;
            dim = std::max(dim, 0);
            dimension = dim;

            simplex_tree = new Simplex_tree_owner(n, size);

            if (DEBUG_FLAG_0) {
                std::cout << points.size() << " vs cntr " << counter << std::endl;

                for (const auto &e: points) {
                    std::cout << e.size() << ": ";
                    for (const auto &a: e) {
                        std::cout << a << ' ';
                    }
                    std::cout << std::endl;
                }
            }
            distances = std::vector<std::vector<double>>(size);
            for (int i = 0; i < distances.size(); ++i) {
                distances[i] = std::vector<double>(i);
                for (int j = 0; j < i; ++j) {
                    double sum_of_squared_differences = 0;
                    for (int k = 0; k < dimension; ++k) {
                        sum_of_squared_differences += (points[i][k] - points[j][k]) * (points[i][k] - points[j][k]);
                    }
                    distances[i][j] = sqrt(sum_of_squared_differences);
                }
            }
        }

        ~Cloud() {
            if (DEBUG_FLAG_0) {
                std::cout << "Destructor Cloud" << std::endl;
            }
            if (simplex_tree == nullptr) {
                std::cout << "NULL" << std::endl;
            }
            simplex_tree->free_tree();
            delete simplex_tree;
        }

        friend std::ostream &operator<<(std::ostream &os, const std::vector<std::vector<double>> &matrix);

    private:
        Cloud(const Cloud &other_cloud) {

        }
    };


}
#endif //CLOUD_H
