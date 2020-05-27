#pragma once
#ifndef SIMPLEX_TREE_OWNER_H
#define SIMPLEX_TREE_OWNER_H

#include "Simplex_tree_node.h"
#include "Simplex_tree_node_inner.h"

namespace smpl {

    struct Simplex_tree_owner {

    private:
        Simplex_tree_node_inner *all_first_vertices;
    public:
        bool to_delete = true;
        std::vector<std::vector<Simplex_tree_node *>> simplices;
    public:

        void free_tree() {
            delete all_first_vertices;
        }

        Simplex_tree_node *get_node_by_position(int pos) {
            int p = 0;
            for (p = 0; p < simplices.size(); ++p) {
                if (simplices[p].size() < pos) {
                    pos -= simplices[p].size();
                } else {
                    break;
                }
            }
            if (p >= simplices.size()) {
                return nullptr;
            }
            return simplices[p][pos];

        }

        Simplex_tree_owner(int n, int cloud_size) {
            all_first_vertices = new Simplex_tree_node_inner();
            simplices = std::vector<std::vector<Simplex_tree_node * >>(n);
            for (int i = 0; i < cloud_size; ++i) {
                auto new_node = new Simplex_tree_node_inner(i, epsilon_a);
                (*(((Simplex_tree_node_inner *) all_first_vertices)->get_successors()))[i] = new_node;
                simplices[0].emplace_back(new_node);
            }
        }

        Simplex_tree_owner(Simplex_tree_owner *simplex_tree, int n) {
            to_delete = false;
            all_first_vertices = simplex_tree->all_first_vertices;
            simplices = std::vector<std::vector<Simplex_tree_node * >>(n);
        }


        Simplex_tree_node *find(const std::vector<int> &simplex, int not_included) {
            Simplex_tree_node *current_simplex_tree_node = all_first_vertices;
            if (current_simplex_tree_node == nullptr) {
                //            unlock_shared();
                return nullptr;
            }
            //        current_simplex_tree_node->lock_shared();
            for (int i = 0; i < simplex.size(); ++i) {
                if (i == not_included) {
                    continue;
                }
                if (i != 0) {
                    ((Simplex_tree_node_inner *) current_simplex_tree_node)->lock_shared();
                }

                Simplex_tree_node *next_node = (*(((Simplex_tree_node_inner *) current_simplex_tree_node)->get_successors()))[simplex[i]];

                if (i != 0) {
                    ((Simplex_tree_node_inner *) current_simplex_tree_node)->unlock_shared();
                }
                current_simplex_tree_node = next_node;
            }
            //        current_simplex_tree_node->unlock_shared();
            return current_simplex_tree_node;
        }

        void insert(Simplex_tree_node *current_simplex_tree_node_temp, int i, double max_birth_time, int vertex) {

            bool is_last = (max_birth_time < 0);
            max_birth_time = abs(max_birth_time);

            auto current_simplex_tree_node = (Simplex_tree_node_inner *) current_simplex_tree_node_temp;
            if (!is_last) {

            }
            current_simplex_tree_node->lock_unique();


            if (is_last) {
                auto new_simplex_node = new Simplex_tree_node(vertex, max_birth_time);

                if (i > 0) {
                    new_simplex_node->set_previous(current_simplex_tree_node);
                }

                if ((current_simplex_tree_node->get_successors()->insert(
                        std::make_pair(vertex, new_simplex_node))).second) {
                    simplices[i].emplace_back(new_simplex_node);
                } else {
                    auto found = current_simplex_tree_node->get_successors()->find(vertex);

                    //                    if (found->second->get_vertex_number() > 0)
                    {
                        simplices[i].emplace_back(found->second);
                        //                    } else {
                        //                        ++extra_cntr;
                    }

                    delete new_simplex_node;
                }
            } else {
                auto new_simplex_node = new Simplex_tree_node_inner(vertex, max_birth_time);

                if (i > 0) {
                    new_simplex_node->set_previous(current_simplex_tree_node);
                }

                if ((current_simplex_tree_node->get_successors()->insert(
                        std::make_pair(vertex, (Simplex_tree_node *) new_simplex_node))).second) {
                    simplices[i].emplace_back((Simplex_tree_node *) new_simplex_node);
                } else {
                    auto found = current_simplex_tree_node->get_successors()->find(vertex);

                    //                    if (found->second->get_vertex_number() > 0)
                    {
                        simplices[i].emplace_back(found->second);
                        //                    } else {
                        //                        ++extra_cntr;
                    }

                    delete new_simplex_node;
                }
            }
            current_simplex_tree_node->unlock_unique();
        }


        void insert(const std::vector<int> &simplex, double birth_time) {

            auto node = find(simplex, simplex.size() - 1);
            insert(node, simplex.size() - 1, birth_time, simplex[simplex.size() - 1]);
        }

        void insert_all_simplices(const std::vector<boost::dynamic_bitset<>> &matrix_of_adjacency,
                                  std::unordered_map<int, int> &new_order_of_points,
                                  const std::vector<int> &subcloud_of_points,
                                  int number_of_vertices,
                                  const std::vector<std::vector<double>> *matrix_of_distances) {
            if (number_of_vertices < 2 || matrix_of_adjacency.empty()) {
                return;
            }

            if (number_of_vertices - 2 >= simplices.size()) {
                return;
            }


            boost::dynamic_bitset<> row(matrix_of_adjacency[0].size());

            for (const auto &temp_subsimplex: simplices[number_of_vertices - 2]) {


                row.set();
                int max_vertex = INT_MIN;

                auto temp_copy_of_subsimplex = temp_subsimplex->get_simplex_by_last_node(number_of_vertices - 1);
                for (const auto &vertex: temp_copy_of_subsimplex) {
                    row &= matrix_of_adjacency[new_order_of_points[vertex]];
                    max_vertex = std::max(new_order_of_points[vertex], max_vertex);
                }
                if (max_vertex != INT_MIN) {
                    auto temp = row.find_next(max_vertex);
                    while (temp < row.size()) {
                        auto found_node = find(temp_copy_of_subsimplex, -1);
                        double max_birth_time = found_node->get_birth_time_protected();

                        for (int i = 0; i < temp_copy_of_subsimplex.size(); ++i) {
                            max_birth_time = std::max(max_birth_time,
                                                      (*matrix_of_distances)[std::max(temp_copy_of_subsimplex[i],
                                                                                      subcloud_of_points[temp])]
                                                      [std::min(temp_copy_of_subsimplex[i], subcloud_of_points[temp])]);
                        }
                        if (number_of_vertices == simplices.size()) {
                            insert(found_node, temp_copy_of_subsimplex.size(), -max_birth_time,
                                   subcloud_of_points[temp]);
                        } else {
                            insert(found_node, temp_copy_of_subsimplex.size(), max_birth_time,
                                   subcloud_of_points[temp]);
                        }
                        temp = row.find_next(temp);
                    }
                }
            }


        }

        ~Simplex_tree_owner() {

        }

    private:
        Simplex_tree_owner(const Simplex_tree_owner &other_simplex_tree) {

        }

    };
}
#endif //SIMPLEX_TREE_OWNER_H
