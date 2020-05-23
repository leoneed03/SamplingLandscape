#pragma once
#ifndef SIMPLEX_TREE_NODE_INNER_H
#define SIMPLEX_TREE_NODE_INNER_H

namespace smpl {

    struct Simplex_tree_node_inner : Simplex_tree_node {
    private:
        //    mutable shared_mutex sh_mutex;
        //    mutable std::mutex sh_mutex;
        tbb::concurrent_unordered_map<int, Simplex_tree_node *> successors;
        //    std::unordered_map<int, Simplex_tree_node *> successors;
        //    boost::container::flat_map<int, Simplex_tree_node *> successors;
    public:

        void lock_unique() {
            //        sh_mutex.lock();
        }

        void unlock_unique() {
            //        sh_mutex.unlock();
        }

        void lock_shared() {
            //        sh_mutex.lock_shared();
            //        sh_mutex.lock();
        }

        void unlock_shared() {
            //        sh_mutex.unlock_shared();
            //        sh_mutex.unlock();
        }


        auto *get_successors() {
            return &successors;
        }

        Simplex_tree_node_inner() {
        }

        Simplex_tree_node_inner(int current_vertex_number, double current_birth_time) {
            vertex_number = current_vertex_number;
            birth_time = current_birth_time;
        }

        ~Simplex_tree_node_inner() {
            for (const auto &e: successors) {
                if (e.second) {
                    if (e.second->get_signed_birth_time() < 0 && abs(e.second->get_birth_time()) > 2 * epsilon_a) {
                        delete (e.second);
                    } else {
                        delete ((Simplex_tree_node_inner *) e.second);
                    }
                }
            }
        }

    private:
        Simplex_tree_node_inner(const Simplex_tree_node_inner &other_simplex_tree_node_inner) {

        }
    };
}

#endif //SIMPLEX_TREE_NODE_INNER_H
