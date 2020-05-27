#pragma once
#ifndef SIMPLEX_TREE_NODE_H
#define SIMPLEX_TREE_NODE_H

#define space_alg 5
#define max_rad 100
#define DEBUG_FLAG_0 (false)
#define epsilon_a (1e-10)

namespace smpl {

//    std::atomic_int zero_cntr = {0};
//    std::atomic_int matrix_size_cntr = {0};
//    std::atomic_int extra_cntr = {0};
//    std::mutex mute;

    struct Simplex_tree_node {
    protected:
        std::atomic_int vertex_number;
        double birth_time = -epsilon_a;
        Simplex_tree_node *previous_vertex = nullptr;

    public:
        double get_signed_birth_time() const {
            return birth_time;
        }

        void set_vertex_number_deleted() {
            vertex_number = -abs(vertex_number);
        }

        int get_abs_vertex_number() const {
            return vertex_number;
        }

        int get_vertex_number() const {
            return vertex_number;
        }

        void set_previous(Simplex_tree_node *simplex_tree_node) {
            previous_vertex = simplex_tree_node;
        }

        double get_birth_time_protected() const {
            return (abs(birth_time) < 2 * epsilon_a) ? (0) : abs(birth_time);
        }

        double get_birth_time() const {
            return (abs(birth_time) < 2 * epsilon_a) ? (0) : abs(birth_time);
        }


        void set_prev(Simplex_tree_node *prev_node) {
            previous_vertex = prev_node;
        }

        std::vector<int> get_simplex_by_last_node(int size) {
            auto current_node = this;
            std::vector<int> result(size);
            int p = size - 1;

            while (current_node != nullptr) {
                result[p] = abs(current_node->vertex_number);
                --p;
                current_node = current_node->previous_vertex;

            }
            return result;
        }

        std::vector<int> get_simplex_by_last_node() {
            auto current_node = this;
            std::vector<int> result;
            while (current_node != nullptr) {
                result.push_back(abs(current_node->vertex_number));
                current_node = current_node->previous_vertex;
            }
            reverse(result.begin(), result.end());
            return result;
        }

        Simplex_tree_node() {

        }

        Simplex_tree_node(int current_vertex_number, double current_birth_time) :
                vertex_number(current_vertex_number),
                birth_time(-current_birth_time) {
        }

    private:
        Simplex_tree_node(const Simplex_tree_node &other_simplex_tree_node) {

        }

    };
}

#endif //SIMPLEX_TREE_NODE_H
