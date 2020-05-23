#pragma once
#ifndef SUBCLOUD_H
#define SUBCLOUD_H

#include "Cloud.h"
#include "Betti_matrix.h"

namespace smpl {
    struct SubCloud {
        std::unordered_map<int, int> new_order_of_points;
        std::vector <boost::dynamic_bitset<>> adjacency_matrix;
        std::vector<int> subcloud_of_points;
        Simplex_tree_owner *root;
        const std::vector <std::vector<double>> *matrix_of_distances;
        double max_radii = 1e6;
        int simplex_counter = 0;
        std::vector <std::vector<int>> new_matrix;
        std::vector<double> from_index_to_birth;
        std::vector<int> dimensions;
        std::unordered_map<Simplex_tree_node *, int> simplex_map;
        std::unordered_map<Simplex_tree_node *, int> filtration_simplex_map;
        std::unordered_map<Simplex_tree_node *, int> less_simplex_map;
        std::unordered_map<Simplex_tree_node *, int> bigger_simplex_map;
        std::vector <std::pair<double, int>> v_pairs;

        ~SubCloud() {
            delete root;
        }

        phat::boundary_matrix <phat::vector_vector> get_boundary_matrix_compressed() {
            phat::boundary_matrix <phat::vector_vector> boundary_matrix;
            std::vector<int> number_of_simplices;
            for (int i = 0; i < root->simplices.size(); ++i) {
                auto bigger_simplices = &root->simplices[i];
                stable_sort(bigger_simplices->begin(), bigger_simplices->end(), [](const auto &lhs, const auto &rhs) {
                    return lhs->get_birth_time() < rhs->get_birth_time();
                });
            }
            int total_size = 0;
            for (int i = 0; i < root->simplices.size(); ++i) {
                total_size += root->simplices[i].size();

                number_of_simplices.push_back(root->simplices[i].size());
            }

            if (DEBUG_FLAG_0) {
                mute.lock();
                std::cout << "simplices: ";
                for (const auto &e: number_of_simplices) {
                    //                fout << e << ' ';
                    std::cout << e << ' ';
                }
                //            fout << '\n';
                std::cout << std::endl << total_size << std::endl;
                mute.unlock();
            }

            boundary_matrix.set_num_cols(total_size);
            v_pairs = std::vector < std::pair < double, int >> (total_size);
            from_index_to_birth = std::vector<double>(total_size);
            dimensions = std::vector<int>(total_size);
            int counter = 0;

            for (int i = 0; i < root->simplices.size(); ++i) {

                auto bigger_simplices = root->simplices[i];
                for (const auto &e: bigger_simplices) {

                    filtration_simplex_map[e] = counter;
                    ++counter;

                    from_index_to_birth[filtration_simplex_map[e]] = e->get_birth_time();
                    dimensions[filtration_simplex_map[e]] = i;
                    v_pairs[filtration_simplex_map[e]] = {e->get_birth_time(), i};

                    boundary_matrix.set_dim(filtration_simplex_map[e], i);
                }

            }
            for (const auto &bigger_simplex_pair: filtration_simplex_map) {
                auto current_bigger_simplex = bigger_simplex_pair.first;
                std::vector <phat::index> temp_col;
                auto current_bigger_simplex_second = current_bigger_simplex->get_simplex_by_last_node();
                int row_number = filtration_simplex_map[current_bigger_simplex];

                if (current_bigger_simplex_second.size() == 1) {
                    continue;
                }

                for (int j = 0; j < current_bigger_simplex_second.size(); ++j) {
                    temp_col.push_back(filtration_simplex_map[root->find(current_bigger_simplex_second, j)]);
                }

                sort(temp_col.begin(), temp_col.end());
                boundary_matrix.set_col(row_number, temp_col);
            }

            return boundary_matrix;
        }

        tbb::concurrent_vector <std::vector<std::pair < double, double>>>

        get_all_dimensions_landscape(bool with_cohomology = false) {

            auto boundary_matrix = get_boundary_matrix_compressed();
            int max_dim = root->simplices.size() - 2;
            if (DEBUG_FLAG_0) {
                std::cout << "Total dim " << max_dim << std::endl;
            }
            tbb::concurrent_vector < std::vector < std::pair < double, double >> > landscape(max_dim + 1);

            phat::chunk_reduction chunk_reduce;
            chunk_reduce(boundary_matrix);
            std::unordered_map<int, int> zero_column_ind;


            std::vector < std::vector < std::pair < int, int >> > p_pairs(max_dim + 1);

            std::vector <phat::persistence_pairs> pairs(max_dim + 1);
            for (int i = 0; i < boundary_matrix.get_num_cols(); ++i) {
                std::vector <phat::index> temp_col;
                boundary_matrix.get_col(i, temp_col);
                if (temp_col.size() == 0 && boundary_matrix.get_dim(i) == max_dim + 1) {
                    ++zero_cntr;
                    auto node = root->get_node_by_position(i);
                    node->set_vertex_number_deleted();
                }
                if (temp_col.size() == 0 && boundary_matrix.get_dim(i) != max_dim + 1) {
                    zero_column_ind[i] = boundary_matrix.get_dim(i);
                    continue;
                }
                if (temp_col.size() == 0) {
                    continue;
                }
                phat::index birth = boundary_matrix.get_max_index(i);
                if (birth < 0) {
                    continue;
                }
                phat::index death = i;
                auto it = zero_column_ind.find(birth);
                if (it != zero_column_ind.end()) {
                    zero_column_ind.erase(it);
                }

                int dimension = boundary_matrix.get_dim(i);
                if (dimension <= 0) {
                    exit(147);
                }
                pairs[dimension - 1].append_pair(birth, death);
            }
            for (const auto &idx: zero_column_ind) {
                if (idx.second == -1 || idx.second == max_dim + 1) {
                    continue;
                }
                p_pairs[idx.second].push_back({idx.first, 888888});
            }

            for (auto &e: pairs) {
                e.sort();
            }
            tbb::concurrent_vector < std::vector < std::pair < double, double >> >
                                                                       persistence_pairs_all_dimensions(max_dim + 1);

            for (int ind = 0; ind < pairs.size(); ++ind) {
                auto &pairr = pairs[ind];
                for (phat::index idx = 0; idx < pairr.get_num_pairs(); idx++) {
                    if (pairr.get_pair(idx).second >= v_pairs.size() || pairr.get_pair(idx).second < 0) {
                        exit(-1);
                    }
                    if (pairr.get_pair(idx).first >= v_pairs.size() || pairr.get_pair(idx).first < 0) {
                        exit(-1);
                    }
                    if (abs(v_pairs[pairr.get_pair(idx).first].first - v_pairs[pairr.get_pair(idx).second].first) >
                        epsilon_a) {
                        persistence_pairs_all_dimensions[ind].push_back({v_pairs[pairr.get_pair(idx).first].first,
                                                                         v_pairs[pairr.get_pair(
                                                                                 idx).second].first});
                    }

                }
            }
            int p = 0;

            for (int p = 0; p < p_pairs.size(); ++p) {
                for (const auto &a: p_pairs[p]) {
                    persistence_pairs_all_dimensions[p].push_back(
                            {v_pairs[a.first].first, std::numeric_limits<double>::max()});
                }
            }

            return persistence_pairs_all_dimensions;
        }

        void get_boundary_matrix_in_column_form() {
            for (int i = 0; i < root->simplices.size(); ++i) {
                auto bigger_simplices = &root->simplices[i];
                stable_sort(bigger_simplices->begin(), bigger_simplices->end(), [](const auto &lhs, const auto &rhs) {
                    return lhs->get_birth_time() < rhs->get_birth_time();
                });
            }
            int total_size = 0;
            for (int i = 0; i < root->simplices.size(); ++i) {
                total_size += root->simplices[i].size();
            }
            from_index_to_birth = std::vector<double>(total_size);
            dimensions = std::vector<int>(total_size);


            for (int i = 0; i < root->simplices.size(); ++i) {

                auto bigger_simplices = root->simplices[i];
                for (const auto &e: bigger_simplices) {
                    simplex_map[e] = simplex_map.size();
                    from_index_to_birth[simplex_map[e]] = e->get_birth_time();
                    dimensions[simplex_map[e]] = i;
                }

            }

            new_matrix = std::vector < std::vector < int >> (simplex_map.size());
            for (const auto &bigger_simplex_pair: simplex_map) {
                auto current_bigger_simplex = bigger_simplex_pair.first;
                auto current_bigger_simplex_second = current_bigger_simplex->get_simplex_by_last_node();
                int row_number = simplex_map[current_bigger_simplex];
                new_matrix[row_number].reserve(current_bigger_simplex_second.size());
                if (current_bigger_simplex_second.size() == 1) {
                    continue;
                }
                for (int j = 0; j < current_bigger_simplex_second.size(); ++j) {
                    new_matrix[row_number].emplace_back(simplex_map[root->find(current_bigger_simplex_second, j)]);
                }
                sort(new_matrix[row_number].begin(), new_matrix[row_number].end());
            }

        }

        void initialize_simplex_tree(int n) {

        }

        void insert_all_simplices_including(int number_of_vertices) {

            for (int i = 3; i <= number_of_vertices; ++i) {
                root->insert_all_simplices(adjacency_matrix, new_order_of_points, subcloud_of_points, i,
                                           matrix_of_distances);
            }
        }

        std::vector <std::vector<std::pair < double, double>>>

        calculate_betti_matrix(int max_number_of_vertices_in_simplex) {
            bool f = true;
            std::vector < std::vector < std::pair < double, double >> > persistence_diagram;

            boost::dynamic_bitset<> already_paired_bigger_simplices(
                    root->simplices[max_number_of_vertices_in_simplex].size(), 0);
            for (int i = 0; i < root->simplices.size(); ++i) {
                auto bigger_simplices = &root->simplices[i];
                stable_sort(bigger_simplices->begin(), bigger_simplices->end(), [](const auto &lhs, const auto &rhs) {
                    return lhs->get_birth_time() < rhs->get_birth_time();
                });
            }

            for (int i = max_number_of_vertices_in_simplex; i >= 1; --i) {


                auto betti_matrix = Betti_matrix(i, &root->simplices[i - 1], &root->simplices[i], root,
                                                 &subcloud_of_points, less_simplex_map, bigger_simplex_map);


                auto persistence_pairs = betti_matrix.construct_betti_matrix(already_paired_bigger_simplices,
                                                                             i == max_number_of_vertices_in_simplex, i);
                std::vector <std::pair<double, double>> persistence_intervals;
                for (const auto &persistence_pair: persistence_pairs) {

                    double birth_time = 0;
                    double death_time;
                    if (persistence_pair.first == INT_MAX) {
                        Simplex_tree_node *stn = ((root->simplices[i])[persistence_pair.second]);
                        birth_time = stn->get_birth_time_protected();
                        persistence_diagram[persistence_diagram.size() - 1].emplace_back(
                                std::make_pair(birth_time, max_radii));
                        death_time = birth_time;
                    } else {
                        if (persistence_pair.second > root->simplices[i - 1].size() - 1) {
                            continue;
                        }
                        auto found_simplex = (root->simplices[i - 1])[root->simplices[i - 1].size() - 1 -
                                                                      persistence_pair.second];
                        auto node_with_this_simplex = (((root->simplices[i - 1])[root->simplices[i - 1].size() - 1 -
                                                                                 persistence_pair.second]));
                        birth_time = node_with_this_simplex->get_birth_time_protected();
                    }
                    if (persistence_pair.first == INT_MAX) {}
                    else {
                        if (persistence_pair.first < 0 || persistence_pair.first >= (root->simplices[i]).size()) {

                        } else {
                            death_time = (((root->simplices[i])[persistence_pair.first]))->get_birth_time_protected();
                        }
                    }
                    if (abs(death_time - birth_time) > epsilon_a) {
                        persistence_intervals.emplace_back(std::make_pair(birth_time, death_time));
                    } else {

                    }
                }

                if (i == 1) {
                    for (int ip = 0; ip < already_paired_bigger_simplices.size(); ++ip) {
                        if (!already_paired_bigger_simplices[ip]) {
                            persistence_intervals.emplace_back(std::make_pair(
                                    ((root->simplices[i - 1])[root->simplices[i - 1].size() - 1 -
                                                              ip])->get_birth_time_protected(), max_radii));
                        }
                    }

                }
                persistence_diagram.emplace_back(persistence_intervals);

            }
            return persistence_diagram;
        }

        SubCloud(const Cloud &cloud, const std::vector<int> &points, double max_radius,
                 int max_number_of_points_in_simplex) {
            max_radii = max_radius;
            matrix_of_distances = &cloud.distances;
            subcloud_of_points = points;
            root = new Simplex_tree_owner(cloud.simplex_tree, cloud.simplex_tree->simplices.size());

            for (const auto &point: points) {
                auto found = root->find(std::vector<int>({point}), -1);
                root->simplices[0].emplace_back(found);
            }


            adjacency_matrix = std::vector<boost::dynamic_bitset<>>(points.size(),
                                                                    boost::dynamic_bitset<>(points.size())); //bitset
            for (const auto &point: points) {
                new_order_of_points[point] = static_cast<int>(new_order_of_points.size());
            }
            initialize_simplex_tree(max_number_of_points_in_simplex);

            for (int i = 0; i < adjacency_matrix.size(); ++i) {
                for (int j = 0; j < adjacency_matrix[i].size(); ++j) {
                    if (points[i] == points[j]) {
                        adjacency_matrix[i][j] = false;
                        continue;
                    }
                    if (cloud.distances[std::max(points[i], points[j])][std::min(points[i], points[j])] < max_radius) {
                        adjacency_matrix[i][j] = true;
                        if (points[i] < points[j]) {
                            root->insert(std::vector<int>({points[i], points[j]}),
                                         cloud.distances[std::max(points[i], points[j])][std::min(points[i],
                                                                                                  points[j])]);
                        }
                    }
                }
            }
        }

    };
}
#endif //SUBCLOUD_H
