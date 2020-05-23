#pragma once
#ifndef BETTI_MATRIX_H
#define BETTI_MATRIX_H

#include "Simplex_tree_owner.h"


namespace smpl {

    void xor_vector(std::vector<int> &l, std::vector<int> &r) {
        std::vector<int> v(l.size() + r.size());
        std::vector<int>::iterator it;

        it = set_symmetric_difference(l.begin(), l.end(), r.begin(), r.end(), v.begin());
        v.resize(it - v.begin());
        auto v1(v);
        swap(l, v1);
    }

    struct Betti_matrix {
        int number_of_points_in_the_less_simplex;
        Simplex_tree_owner *simplex_tree;
        std::vector<Simplex_tree_node *> *less_simplices;
        std::vector<Simplex_tree_node *> *bigger_simplices;
        std::vector<int> *original_numbers_of_points;
        std::unordered_map<Simplex_tree_node *, int> &less_simplex_map_;

        Betti_matrix(int new_number_of_points_in_the_less_simplex,
                     std::vector<Simplex_tree_node *> *new_less_simplices,
                     std::vector<Simplex_tree_node *> *new_bigger_simplices,
                     Simplex_tree_owner *new_simplex_tree,
                     std::vector<int> *new_original_numbers_of_points,
                     std::unordered_map<Simplex_tree_node *, int> &new_less_simplex_map,
                     std::unordered_map<Simplex_tree_node *, int> &new_bigger_simplex_map) :
                simplex_tree(new_simplex_tree),
                original_numbers_of_points(new_original_numbers_of_points),
                number_of_points_in_the_less_simplex(new_number_of_points_in_the_less_simplex),
                less_simplices(new_less_simplices),
                bigger_simplices(new_bigger_simplices),
                less_simplex_map_(new_less_simplex_map) {}

        //            rows are written from right to left
        //            rows are bigger simplices
        //            columns are less
        //            bits is "should be used" are inversed (columns are from lowest to biggest number and row vice versa)!!!!!!!!!!!!!!!!!!!!!!!!
        //
        std::vector <std::vector<int>> *new_matrix_pointer = nullptr;


        std::vector <std::pair<int, int>> reduce_matrix(std::vector <std::vector<int>> &temp_betti_matrix,
                                                        boost::dynamic_bitset<> &paired_bigger_simplices_should_not_be_used,
                                                        bool is_max_simplex_rank,
                                                        std::priority_queue<int, std::vector<int>, std::function<bool(
                                                                const int &,
                                                                const int &)>> &rows_with_lowest_bits) {
            int xor_counter = 0;
            int max_size = 0;
            int mm_size = 0;
            for (const auto &e: temp_betti_matrix) {
                max_size += e.size();
            }


            std::vector <std::pair<int, int>> persistence_pairs;
            int max_number_of_vertices1 = INT_MIN, max_number_of_vertices2 = INT_MIN;
            persistence_pairs.reserve(std::min(less_simplices->size(), bigger_simplices->size()));
            int pos_in_paired_simplices = paired_bigger_simplices_should_not_be_used.find_first();    //this is reversed std::vector
            boost::dynamic_bitset<> should_use(temp_betti_matrix.size());
            should_use.set(); //while normalizing matrix we have to mark used rows in order not to normalize them again
            if (!is_max_simplex_rank) {
                while (pos_in_paired_simplices < paired_bigger_simplices_should_not_be_used.size()) {
                    ////            betti_matrix[pos_in_paired_simplices].reset();
                    temp_betti_matrix[pos_in_paired_simplices] = {};
                    should_use[pos_in_paired_simplices] = 0;
                    pos_in_paired_simplices = paired_bigger_simplices_should_not_be_used.find_next(
                            pos_in_paired_simplices);
                }

            }
            boost::dynamic_bitset<> temp_paired_simplices;
            if (!is_max_simplex_rank) {
                temp_paired_simplices = boost::dynamic_bitset<>(paired_bigger_simplices_should_not_be_used);
            }
            paired_bigger_simplices_should_not_be_used = boost::dynamic_bitset<>(less_simplices->size(),
                                                                                 0); //we will pass this bitset consisting of already paired simplices to next iteration
            if (temp_betti_matrix.size() == 0) {
                return persistence_pairs;
            }

            int lowest_row = INT_MAX;


            int max_counter = INT_MIN;
            while (!rows_with_lowest_bits.empty()) {

                lowest_row = rows_with_lowest_bits.top();
                int lowest_bit = temp_betti_matrix[lowest_row][0];
                rows_with_lowest_bits.pop();
                if (should_use[lowest_row]) {
                    should_use[lowest_row] = 0;
                    persistence_pairs.emplace_back(std::make_pair(lowest_row, lowest_bit));
                    paired_bigger_simplices_should_not_be_used[paired_bigger_simplices_should_not_be_used.size() - 1 -
                                                               lowest_bit] = 1;
                }
                if (rows_with_lowest_bits.empty()) {
                    break;
                }
                auto next_top_bit_and_row = rows_with_lowest_bits.top();
                while (!rows_with_lowest_bits.empty() &&
                       temp_betti_matrix[rows_with_lowest_bits.top()][0]/*rows_with_lowest_bits.top().first*/ ==
                       lowest_bit) {
                    next_top_bit_and_row = rows_with_lowest_bits.top();
                    rows_with_lowest_bits.pop();

                    max_size -= temp_betti_matrix[next_top_bit_and_row].size();
                    max_size -= temp_betti_matrix[lowest_row].size();
                    auto temp_max_counter = std::max(temp_betti_matrix[next_top_bit_and_row].size(),
                                                     temp_betti_matrix[lowest_row].size());

                    xor_vector(temp_betti_matrix[next_top_bit_and_row], temp_betti_matrix[lowest_row]);
                    auto temp2_max_counter = std::max(temp_betti_matrix[next_top_bit_and_row].size(),
                                                      temp_betti_matrix[next_top_bit_and_row].size());

                    max_size += temp_betti_matrix[next_top_bit_and_row].size();
                    max_size += temp_betti_matrix[lowest_row].size();
                    mm_size = std::max(mm_size, max_size);
                    max_counter = std::max(max_counter, (int) std::max(temp2_max_counter, temp_max_counter));

                    ++xor_counter;
                    if (temp_betti_matrix[next_top_bit_and_row].size() > 0) {
                        rows_with_lowest_bits.emplace(next_top_bit_and_row);
                    }

                }

                std::vector<int> t(1, lowest_bit);
                swap(t, temp_betti_matrix[lowest_row]);

            }

            if (!is_max_simplex_rank) {
                for (int i = 0; i < temp_betti_matrix.size(); ++i) {
                    if (!is_max_simplex_rank && temp_betti_matrix[i].empty() && !temp_paired_simplices[i]) {
                        persistence_pairs.emplace_back(std::make_pair(INT_MAX, i));
                    }
                }

            }
            return persistence_pairs;
        }

        std::vector <std::pair<int, int>> construct_betti_matrix(boost::dynamic_bitset<> &already_paired_simplices,
                                                                 bool is_max_simplex_rank,
                                                                 int current_number_of_vertices) {
            bool f = true;
            int simplex_indicator;

            simplex_indicator = 0;
            std::unordered_map < Simplex_tree_node * , int > less_simplex_map;
            less_simplex_map.reserve(less_simplices->size());
            for (const auto &current_less_simplex: *less_simplices) {
                less_simplex_map.emplace(std::make_pair(current_less_simplex, simplex_indicator));
                ++simplex_indicator;
            }

            std::vector <std::vector<int>> new_matrix(bigger_simplices->size());

            std::vector<int> vector_for_pq;

            if (is_max_simplex_rank) {
                vector_for_pq.reserve(bigger_simplices->size());
            }
            new_matrix_pointer = &new_matrix;
            for (int i = 0; i < bigger_simplices->size(); ++i) {
                auto &current_bigger_simplex = (*bigger_simplices)[i];
                auto current_bigger_simplex_second = current_bigger_simplex->get_simplex_by_last_node(
                        current_number_of_vertices + 1);
                int row_number = i;
                new_matrix[row_number].reserve(current_bigger_simplex_second.size());
                for (int j = 0; j < current_bigger_simplex_second.size(); ++j) {
                    new_matrix[row_number].emplace_back(less_simplices->size() - 1 -
                                                        less_simplex_map[simplex_tree->find(
                                                                current_bigger_simplex_second,
                                                                j)]);
                }
                sort(new_matrix[row_number].begin(), new_matrix[row_number].end());
                vector_for_pq.emplace_back(row_number);
            }


            std::priority_queue<int, std::vector<int>, std::function<bool(const int &,
                                                                          const int &)>> rows_with_lowest_bits{
                    [&](const int &rhs, const int &lhs) -> bool {
                        return (*new_matrix_pointer)[lhs][0] < (*new_matrix_pointer)[rhs][0] ||
                               (!((*new_matrix_pointer)[rhs][0] < (*new_matrix_pointer)[lhs][0]) && lhs < rhs);
                    }, std::move(vector_for_pq)};

            return reduce_matrix(new_matrix, already_paired_simplices, is_max_simplex_rank, rows_with_lowest_bits);
        }

    };
}
#endif //BETTI_MATRIX_H
