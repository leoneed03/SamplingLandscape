#include <phat/compute_persistence_pairs.h>
#include <phat/algorithms/chunk_reduction.h>

#include <atomic>
#include <iostream>
#include <vector>
#include <algorithm>
#include <set>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <map>
#include <thread>
#include <shared_mutex>
#include <mutex>
#include <queue>

#include <tbb/concurrent_unordered_map.h>
#include <tbb/concurrent_vector.h>
#include <boost/container/flat_map.hpp>
#include <unordered_map>
#include <chrono>
#include <random>

#include <gudhi/Persistence_landscape.h>

/*
#include <phat/compute_persistence_pairs.h>
#include <phat/representations/vector_vector.h>
#include <phat/algorithms/chunk_reduction.h>

*/
//#include <phat/helpers/dualize.h>
//#include <phat/representations/heap_pivot_column.h>
//#include <phat/algorithms/standard_reductistrion.h>
//#include <phat/algorithms/row_reduction.h>
//#include <phat/algorithms/twist_reduction.h>


//#include "../../../phat/compute_persistence_pairs.h"
//#include "../../../phat/representations/vector_vector.h"
//#include "../../../phat/algorithms/chunk_reduction.h"

#include <boost/bind.hpp>
#include <boost/thread/thread.hpp>
#include <boost/asio.hpp>
#include <boost/dynamic_bitset.hpp>

#include "../mean_landscapes/mean_landscapes.cpp"
#define space 5
#define max_rad 100

#define epsilon_a (1e-10)

//using namespace std;


std::atomic_int zero_cntr = std::atomic_int(0);
std::atomic_int matrix_size_cntr = std::atomic_int(0);
std::atomic_int extra_cntr = std::atomic_int(0);
std::mutex mute;


bool flag = false;


int number_of_dots_in_file;
//int number_of_errors = 0;

using Persistence_landscape = Gudhi::Persistence_representations::Persistence_landscape;

std::ostream &operator<<(std::ostream &os, const std::vector <boost::dynamic_bitset<>> &matrix) {
    os << std::endl;
    for (int j = 0; j < matrix.size(); ++j) {
        boost::dynamic_bitset<> row = matrix[j];
        os << std::setw(space) << j << ": ";
        for (int i = 0; i < row.size(); ++i) {
            os << row[i];
        }
        os << std::endl;
    }
    return os;
}

void xor_vector(std::vector<int> &l, std::vector<int> &r) {
    std::vector<int> v(l.size() + r.size());                      // 0  0  0  0  0  0  0  0  0  0
    std::vector<int>::iterator it;

    it = set_symmetric_difference(l.begin(), l.end(), r.begin(), r.end(), v.begin());
    v.resize(it - v.begin());
//    std::vector<int>
    auto v1(v);
//    l = std::move(v);
//    auto vector_copy(std::move(v));
    swap(l, v1);
}


struct Simplex_tree_node {
protected:
    std::atomic_int vertex_number;
    double birth_time = -epsilon_a;
    Simplex_tree_node *previous_vertex = nullptr;

public:
    double get_signed_birth_time() {
        return birth_time;
    }

    void set_vertex_number_deleted() {
        vertex_number = -abs(vertex_number);
    }

    int get_abs_vertex_number() {
        return vertex_number;
    }

    int get_vertex_number() {
        return vertex_number;
    }

    void set_previous(Simplex_tree_node *simplex_tree_node) {
        previous_vertex = simplex_tree_node;
    }

    double get_birth_time_protected() {
        return (abs(birth_time) < 2 * epsilon_a) ? (0) : abs(birth_time);
    }

    double get_birth_time() {
        return (abs(birth_time) < 2 * epsilon_a) ? (0) : abs(birth_time);
    }


    void set_prev(Simplex_tree_node *prev_node) {
        previous_vertex = prev_node;
    }


//    std::vector<int> get_simplex_by_last_node() {
//        auto current_node = this;
//        std::vector<int> result;
//        while (current_node != nullptr) {
//            result.emplace_back(current_node->vertex_number);
////
//            current_node->previous_vertex->sh_mutex.lock_shared();
//            auto prev = current_node;
//            current_node = current_node->previous_vertex;
//            prev->sh_mutex.unlock_shared();
//        }
//        reverse(result.begin(), result.end());
//        return result;
//    }

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

};


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

    /*
    void print(int level) {
        for (int i = 0; i < level; ++i) {
             //std::cout << "*";
        }
         //std::cout << " " << get_vertex_number() << " [" << this->get_birth_time() << "]" << std::endl;
        for (const auto &e: successors) {
            if (e.second->get_signed_birth_time() < 0) {
                for (int i = 0; i <= level; ++i) {
                     //std::cout << "*";
                }
                 //std::cout << " " << e.second->get_vertex_number() << " [" << e.second->get_signed_birth_time() << "]"
                     << std::endl;
            } else {
                ((Simplex_tree_node_inner *) e.second)->print(level + 1);
            }
        }
    }*/
    /*
    void print() {
         //std::cout << std::endl;
        print(1);
    }*/
};

struct Simplex_tree {     //this is only 1's thread property

public:
    Simplex_tree_node_inner *all_first_vertices; //shared resource
public:
    bool to_delete = true;
    std::vector <std::vector<Simplex_tree_node *>> simplices; //private for thread-owner


public:
//    void lock_unique() {
//        sh_mutex.lock();
//    }
//    void unlock_unique() {
//        sh_mutex.unlock();
//    }
//    void lock_shared() {
//        sh_mutex.lock_shared();
//    }
//    void unlock_shared() {
//        sh_mutex.unlock_shared();
//    }

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

    Simplex_tree(int n, int cloud_size) {
        all_first_vertices = new Simplex_tree_node_inner();
        simplices = std::vector < std::vector < Simplex_tree_node * >> (n);
        for (int i = 0; i < cloud_size; ++i) {
            auto new_node = new Simplex_tree_node_inner(i, epsilon_a);
            (*(((Simplex_tree_node_inner *) all_first_vertices)->get_successors()))[i] = new_node;
            simplices[0].emplace_back(new_node);
        }
    }

    Simplex_tree(Simplex_tree *simplex_tree, int n) {
        to_delete = false;
        all_first_vertices = simplex_tree->all_first_vertices;
        simplices = std::vector < std::vector < Simplex_tree_node * >> (n);
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
            //changed
            current_simplex_tree_node = next_node;

//            if (current_simplex_tree_node == nullptr) {
//                return nullptr;
//            }
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

            if ((current_simplex_tree_node->get_successors()->insert(std::make_pair(vertex, new_simplex_node))).second) {
                simplices[i].emplace_back(new_simplex_node);
            } else {
                auto found = current_simplex_tree_node->get_successors()->find(vertex);

                if (found->second->get_vertex_number() > 0) {
                    simplices[i].emplace_back(found->second);
                } else {
                    ++extra_cntr;
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
                if (found->second->get_vertex_number() > 0) {
                    simplices[i].emplace_back(found->second);  //here

                } else {
                    ++extra_cntr;
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

    void insert_all_simplices(const std::vector <boost::dynamic_bitset<>> &matrix_of_adjacency,
                              std::unordered_map<int, int> &new_order_of_points,
                              const std::vector<int> &subcloud_of_points,
                              int number_of_vertices,
                              const std::vector <std::vector<double>> *matrix_of_distances) {
        if (number_of_vertices < 2 || matrix_of_adjacency.empty()) {
            return;
        }
        if (number_of_vertices - 2 >= simplices.size()) {
            return;
        }

//        for (const auto &e: subcloud_of_points) {
//
//        }


        boost::dynamic_bitset<> row(matrix_of_adjacency[0].size());

        for (const auto &temp_subsimplex: simplices[number_of_vertices - 2]) {


            row.set();
            int max_vertex = INT_MIN;

            auto temp_copy_of_subsimplex = temp_subsimplex->get_simplex_by_last_node(number_of_vertices - 1);
//            for (const auto& e: temp_copy_of_subsimplex) {
//
//            }
//
//
//            if (temp_subsimplex2.second.size() > number_of_vertices - 1)
//                exit(17);
            for (const auto &vertex: temp_copy_of_subsimplex) {
                row &= matrix_of_adjacency[new_order_of_points[vertex]];
                max_vertex = std::max(new_order_of_points[vertex], max_vertex);
            }
            if (max_vertex != INT_MIN) {
                auto temp = row.find_next(max_vertex);
                while (temp < row.size()) {
//                    temp_copy_of_subsimplex = temp_subsimplex2.second;
                    auto found_node = find(temp_copy_of_subsimplex, -1);
                    double max_birth_time = found_node->get_birth_time_protected();

                    for (int i = 0; i < temp_copy_of_subsimplex.size(); ++i) {
                        max_birth_time = std::max(max_birth_time,
                                             (*matrix_of_distances)[std::max(temp_copy_of_subsimplex[i],
                                                                        subcloud_of_points[temp])]
                                             [std::min(temp_copy_of_subsimplex[i], subcloud_of_points[temp])]);
                    }
                    if (number_of_vertices == simplices.size()) {
                        insert(found_node, temp_copy_of_subsimplex.size(), -max_birth_time, subcloud_of_points[temp]);
                    } else {
                        insert(found_node, temp_copy_of_subsimplex.size(), max_birth_time, subcloud_of_points[temp]);
                    }
                    temp = row.find_next(temp);
                }
            }
        }


    }

//    void print() {
//        (all_first_vertices)->print();
//    }
    void print() {
        for (const auto &e: *(((Simplex_tree_node_inner *) all_first_vertices)->get_successors())) {

            if (e.second) {
                //((Simplex_tree_node_inner *) e.second)->print();
            }
        }
//        for (const auto &simplices_array: simplices) {
//            for (const auto &simplex: simplices_array) {
//
//                for (const auto &e: simplex.second->get_simplex_by_last_node()) {
//
//                }
//
//            }
//        }
    }

    ~Simplex_tree() {
//        if (to_delete) {
//            delete ((Simplex_tree_node_inner *) all_first_vertices);
//            all_first_vertices = nullptr;
//        }
    }

};


struct Cloud {
    int dimension, size;
    std::vector <std::vector<double>> points;
    std::vector <std::vector<double>> distances;
    Simplex_tree *simplex_tree;

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

//        std::cout << "_" << line << std::endl;
//        std::getline(input_stream, line);
//        std::cout << "_" << line << std::endl;

        while (std::getline(input_stream, line)) {
//            std::cout << "_________" << line << std::endl;
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
            assert(point.size() == points.front().size());
        }






//        fin >> size >> dimension;
        size = counter;
        dimension = dim;
        simplex_tree = new Simplex_tree(n, size);

        std::cout << points.size() << " vs cntr " << counter << std::endl;
        for (const auto& e: points) {
            std::cout << e.size() << ": ";
            for (const auto& a: e) {
                std::cout << a << ' ';
            }
            std::cout << std::endl;
        }
//        simplex_tree->print();
//         //std::cout << "printed|n\n\n\n" << std::endl;
//        number_of_dots_in_file = size;

//        points = std::vector < std::vector < double >> (size, std::vector<double>(dimension, 0));
//        for (int i = 0; i < points.size(); ++i) {
//            for (int j = 0; j < dimension; ++j) {
//                fin >> points[i][j];
//            }
//        }
        distances = std::vector < std::vector < double >> (size);
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
//        for (int i = 0; i < points.size(); ++i) {
//            simplex_tree->insert(std::vector<int>({i}), 0);
//        }
//        for (const auto &row: distances) {
//            for (const auto &e: row) {
//                 //std::cout << std::setw(space + 5) << e;
//            }
//             //std::cout << std::endl;
//        }

    }

    ~Cloud() {
        std::cout << "Destructor Cloud" << std::endl;
        if (simplex_tree == nullptr) {
            std::cout << "NULL" << std::endl;
        }
        delete simplex_tree->all_first_vertices;
        delete simplex_tree;
    }

    friend std::ostream &operator<<(std::ostream &os, const std::vector <std::vector<double>> &matrix);

    void print_distances(std::ostream &os) {

        for (const auto &row: distances) {
            for (const auto &element: row) {

            }

        }
    }
};


struct Betti_matrix {
    int number_of_points_in_the_less_simplex;
    Simplex_tree *simplex_tree;
    std::vector<Simplex_tree_node *> *less_simplices;
    std::vector<Simplex_tree_node *> *bigger_simplices;
    std::vector<int> *original_numbers_of_points;

    std::unordered_map<Simplex_tree_node *, int> &less_simplex_map_;
//    std::unordered_map<Simplex_tree_node *, int>& bigger_simplex_map;

    Betti_matrix(int new_number_of_points_in_the_less_simplex,
                 std::vector<Simplex_tree_node *> *new_less_simplices,
                 std::vector<Simplex_tree_node *> *new_bigger_simplices,
                 Simplex_tree *new_simplex_tree,
                 std::vector<int> *new_original_numbers_of_points,
                 std::unordered_map<Simplex_tree_node *, int> &new_less_simplex_map,
                 std::unordered_map<Simplex_tree_node *, int> &new_bigger_simplex_map) :
            simplex_tree(new_simplex_tree),
            original_numbers_of_points(new_original_numbers_of_points),
            number_of_points_in_the_less_simplex(new_number_of_points_in_the_less_simplex),
            less_simplices(new_less_simplices),
            bigger_simplices(new_bigger_simplices),
            less_simplex_map_(new_less_simplex_map)
    /*bigger_simplex_map(new_bigger_simplex_map)*/ {}

//            rows are written from right to left
//            rows are bigger simplices
//            columns are less
//            bits is "should be used" are inversed (columns are from lowest to biggest number and row vice versa)!!!!!!!!!!!!!!!!!!!!!!!!
//
    std::vector <std::vector<int>> *new_matrix_pointer = nullptr;


    std::vector <std::pair<int, int>> reduce_matrix(std::vector <std::vector<int>> &temp_betti_matrix,
                                          boost::dynamic_bitset<> &paired_bigger_simplices_should_not_be_used,
                                          bool is_max_simplex_rank,
                                          std::priority_queue<int, std::vector<int>, std::function<bool(const int &,
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
//        int pos_in_paired_simplices = paired_bigger_simplices_should_be_used.find_first();
        boost::dynamic_bitset<> should_use(temp_betti_matrix.size());
//        boost::dynamic_bitset<> paired_simplices(betti_matrix.size());
//        paired_simplices.std::set();
        should_use.set(); //while normalizing matrix we have to mark used rows in order not to normalize them again
        if (!is_max_simplex_rank) {
            while (pos_in_paired_simplices < paired_bigger_simplices_should_not_be_used.size()) {
                ////            betti_matrix[pos_in_paired_simplices].reset();
                temp_betti_matrix[pos_in_paired_simplices] = {};
                should_use[pos_in_paired_simplices] = 0;
                pos_in_paired_simplices = paired_bigger_simplices_should_not_be_used.find_next(pos_in_paired_simplices);
            }

        }
        boost::dynamic_bitset<> temp_paired_simplices;
        if (!is_max_simplex_rank) {
            temp_paired_simplices = boost::dynamic_bitset<>(paired_bigger_simplices_should_not_be_used);
        }
        paired_bigger_simplices_should_not_be_used = boost::dynamic_bitset<>(less_simplices->size(),
                                                                             0); //we will pass this bitset consisting of already paired simplices to next iteration
        //but it has to be resized (to size of less simplices)
//        paired_bigger_simplices_should_be_used.std::set();

        if (temp_betti_matrix.size() == 0) {
            return persistence_pairs;
        }
//        std::pair<int, int> lowest_bit_and_row = {INT_MAX, INT_MAX};
        int lowest_row = INT_MAX;


        int max_counter = INT_MIN;
        while (!rows_with_lowest_bits.empty()) {

//            lowest_bit_and_row = rows_with_lowest_bits.top();
            lowest_row = rows_with_lowest_bits.top();
            int lowest_bit = temp_betti_matrix[lowest_row][0];
            rows_with_lowest_bits.pop();
            if (should_use[lowest_row]) {
                should_use[lowest_row] = 0;
//                std::vector<int> t(1, lowest_bit);
//                swap(t, temp_betti_matrix[lowest_row]);
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
//                if (temp_betti_matrix[next_top_bit_and_row.second].size() > 4 &&
//                    (temp_betti_matrix[next_top_bit_and_row.second][3] ==
//                     temp_betti_matrix[next_top_bit_and_row.second][2] ||
//                     temp_betti_matrix[next_top_bit_and_row.second][1] ==
//                     temp_betti_matrix[next_top_bit_and_row.second][2]))
//                    exit(111);
//
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

//        boost::container::flat_map<Simplex_tree_node*, int> less_simplex_map;
        std::unordered_map < Simplex_tree_node * , int > less_simplex_map;
        less_simplex_map.reserve(less_simplices->size());
        for (const auto &current_less_simplex: *less_simplices) {
//            simplex_map.emplace(std::make_pair(current_less_simplex, simplex_indicator));

            less_simplex_map.emplace(std::make_pair(current_less_simplex, simplex_indicator));
            ++simplex_indicator;
        }
//        if (is_max_simplex_rank) {
//            simplex_indicator = 0;
//            bigger_simplex_map.reserve(bigger_simplices->size());
//            for (const auto &current_bigger_simplex: *bigger_simplices) {
//                bigger_simplex_map.emplace(std::make_pair(current_bigger_simplex, simplex_indicator));
//                ++simplex_indicator;
//            }
//        }


        std::vector <std::vector<int>> new_matrix(bigger_simplices->size());


        std::vector<int> vector_for_pq;


        if (is_max_simplex_rank) {
            vector_for_pq.reserve(bigger_simplices->size());
        }
        new_matrix_pointer = &new_matrix;
        for (int i = 0; i < bigger_simplices->size(); ++i) {
//            if (simplex_map.find(current_bigger_simplex.second) == simplex_map.end()) {
//                exit(1);
//            }
            auto &current_bigger_simplex = (*bigger_simplices)[i];
            auto current_bigger_simplex_second = current_bigger_simplex->get_simplex_by_last_node(
                    current_number_of_vertices + 1);
//            int row_number = bigger_simplex_map[current_bigger_simplex];
            int row_number = i;
            new_matrix[row_number].reserve(current_bigger_simplex_second.size());
            for (int j = 0; j < current_bigger_simplex_second.size(); ++j) {
                new_matrix[row_number].emplace_back(less_simplices->size() - 1 -
                                                    less_simplex_map[simplex_tree->find(current_bigger_simplex_second,
                                                                                        j)]);
            }
            sort(new_matrix[row_number].begin(), new_matrix[row_number].end());
            vector_for_pq.emplace_back(row_number);
        }


        std::priority_queue<int, std::vector<int>, std::function<bool(const int &, const int &)>> rows_with_lowest_bits{
                [&](const int &rhs, const int &lhs) -> bool {
//                    return lhs.first > rhs.first;
//                        if (lhs.first > rhs.first) {
//                            return true;
//                        } else {
//                            return lhs.second >= rhs.second;
//                        }
//                    return rhs > lhs;
///////////////smth is wrong - infinite loop
                    return (*new_matrix_pointer)[lhs][0] < (*new_matrix_pointer)[rhs][0] ||
                           (!((*new_matrix_pointer)[rhs][0] < (*new_matrix_pointer)[lhs][0]) && lhs < rhs);
//                        return rhs < lhs;
                }, std::move(vector_for_pq)};

        std::ofstream fout("m" + std::to_string(current_number_of_vertices) + ".sms");
        fout << bigger_simplices->size() << " " << less_simplices->size() << " M\n";
        for (int i = 0; i < new_matrix.size(); ++i) {
            int row = i + 1;
            for (auto &e: new_matrix[i]) {
                fout << row << " " << e + 1 << " 1\n";
            }
        }
        fout << "0 0 0";
        return reduce_matrix(new_matrix, already_paired_simplices, is_max_simplex_rank, rows_with_lowest_bits);
    }

};

struct SubCloud {
    std::unordered_map<int, int> new_order_of_points;
    std::vector <boost::dynamic_bitset<>> adjacency_matrix;
    std::vector<int> subcloud_of_points;
    Simplex_tree *root;
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
//        delete root;
    }

//    template<typename T>
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

        {
            mute.lock();

//            std::ofstream fout("bunny_log.txt");
//            fout << "simplices: ";
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
//        std::set<int> indices;
        int counter = 0;
        for (int i = 0; i < root->simplices.size(); ++i) {

            auto bigger_simplices = root->simplices[i];
            for (const auto &e: bigger_simplices) {

//                if (indices.find(filtration_simplex_map.size()) != indices.end()) {
//
//                     //std::cout << "ooooh" << std::endl;
//                    return boundary_matrix;
//                    exit(14);
//                }
                filtration_simplex_map[e] = counter;
                ++counter;

//                indices.insert(filtration_simplex_map[e]);
                from_index_to_birth[filtration_simplex_map[e]] = e->get_birth_time();
                dimensions[filtration_simplex_map[e]] = i;
                v_pairs[filtration_simplex_map[e]] = {e->get_birth_time(), i};

                boundary_matrix.set_dim(filtration_simplex_map[e], i);
            }

        }
        for (const auto &bigger_simplex_pair: filtration_simplex_map) {
            auto current_bigger_simplex = bigger_simplex_pair.first;
//            for (const auto& e: current_bigger_simplex->get_simplex_by_last_node()) {
//                 //std::cout << e << ' ';
//            }
//             //std::cout << std::endl;
            std::vector <phat::index> temp_col;
            auto current_bigger_simplex_second = current_bigger_simplex->get_simplex_by_last_node();
            int row_number = filtration_simplex_map[current_bigger_simplex];

            if (current_bigger_simplex_second.size() == 1) {
                continue;
            }

            for (int j = 0; j < current_bigger_simplex_second.size(); ++j) {
                temp_col.push_back(filtration_simplex_map[root->find(current_bigger_simplex_second, j)]);
            }
//             //std::cout << " the size is " << temp_col.size() << " for " << row_number << std::endl;


            sort(temp_col.begin(), temp_col.end());
            boundary_matrix.set_col(row_number, temp_col);
        }
//        std::ofstream fout("bm1.txt");
//        bool f = false;
//        for (int i = 0; i < boundary_matrix.get_num_cols(); ++i) {
//            std::vector< phat::index > temp_col;
//            boundary_matrix.get_col(i, temp_col);
//            if (!temp_col.empty()) {
//                f = true;
//            }
//            phat::index birth = boundary_matrix.get_max_index(i);
////            for (const auto& e: temp_col) {
////                fout << e << ' ';
////            }
//            fout << birth << "_\n";
//        }
//        if (!f) {
//            exit(2);
//        }
//         //std::cout << "finished writing" << std::endl;

        return boundary_matrix;
    }


//    template<typename T, typename U>
    tbb::concurrent_vector <std::vector<std::pair < double, double>>>

    get_all_dimensions_landscape(bool with_cohomology = false) {
        auto boundary_matrix = get_boundary_matrix_compressed();

        int max_dim = root->simplices.size() - 2;
        tbb::concurrent_vector < std::vector < std::pair < double, double >> > landscape(max_dim + 1);


//        with_cohomology = false;

        //////block cohomology starts

        if (with_cohomology) {
            phat::persistence_pairs pers_pairs_temp;
            std::cout << "cohomology" << std::endl;
            phat::compute_persistence_pairs_dualized<phat::twist_reduction, phat::vector_vector>(pers_pairs_temp,
                                                                                                 boundary_matrix);
            for (int idx = 0; idx < pers_pairs_temp.get_num_pairs(); ++idx) {
                if (abs(v_pairs[pers_pairs_temp.get_pair(idx).first].first - v_pairs[pers_pairs_temp.get_pair(
                        idx).second].first) > epsilon_a) {
                    landscape[dimensions[pers_pairs_temp.get_pair(idx).first]].push_back(
                            {v_pairs[pers_pairs_temp.get_pair(idx).first].first,
                             v_pairs[pers_pairs_temp.get_pair(
                                     idx).second].first});
                }

            }

            return landscape;
        }
        //////block cohomology ends


        phat::chunk_reduction chunk_reduce;
        chunk_reduce(boundary_matrix);
        std::unordered_map<int, int> zero_column_ind;


        std::vector < std::vector < std::pair < int, int >> > p_pairs(max_dim + 1);

        std::vector <phat::persistence_pairs> pairs(max_dim + 1);
        for (int i = 0; i < boundary_matrix.get_num_cols(); ++i) {
//         //std::cout << i << std::endl;
            std::vector <phat::index> temp_col;
            boundary_matrix.get_col(i, temp_col);
            if (temp_col.size() == 0 && boundary_matrix.get_dim(i) == max_dim + 1) {
//                std::cout << i << " this dimension is " /*<< boundary_matrix.get_dim(i) <<  " with " */<< v_pairs[i].second <<std::endl;
                ++zero_cntr;
                auto node = root->get_node_by_position(i);
                node->set_vertex_number_deleted();
//                if (node->get_vertex_number() < 0) {
//                    std::cout << "RIP";
//                    exit(15);
//                }
//                std::cout << zero_cntr << std::endl;
//                exit(101);
            }
            if (temp_col.size() == 0 && boundary_matrix.get_dim(i) != max_dim + 1) {
//            ++counter;
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
                //std::cout << "dim problem" << std::endl;
                exit(147);
            }
            pairs[dimension - 1].append_pair(birth, death);
//         //std::cout << '[' << birth << ',' << death << ']' << std::endl;
        }
        for (const auto &idx: zero_column_ind) {
            if (idx.second == -1 || idx.second == max_dim + 1) {
                continue;
            }
            p_pairs[idx.second].push_back({idx.first, 888888});
        }


//        for (int i = 0; i < p_pairs.size(); ++i) {
//            for (const auto &e: p_pairs[i]) {
//            }
//        }
        for (auto &e: pairs) {
            e.sort();
        }
        tbb::concurrent_vector < std::vector < std::pair < double, double >> > persistence_pairs_all_dimensions(max_dim + 1);

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
//                    pp.push_back({v_pairs[pairr.get_pair(idx).first].second, {v_pairs[pairr.get_pair(idx).first].first,
//                                                                              v_pairs[pairr.get_pair(
//                                                                                      idx).second].first}});

                    persistence_pairs_all_dimensions[ind].push_back({v_pairs[pairr.get_pair(idx).first].first,
                                                                     v_pairs[pairr.get_pair(
                                                                             idx).second].first});
                }

            }
        }
//        sort(pp.begin(), pp.end());
        int p = 0;


        //std::cout << "do it!" << std::endl;
        for (int p = 0; p < p_pairs.size(); ++p) {
            for (const auto &a: p_pairs[p]) {
//                     //std::cout << e.first << ':' << v_pairs[a.first].first << " inf" << std::endl;
                persistence_pairs_all_dimensions[p].push_back(
                        {v_pairs[a.first].first, std::numeric_limits<double>::max()});
//                if (max_radii > 1e10) {
//                    //std::cout << "dead" << std::endl;
//                    exit(15);
//                }
            }
//                ++p;


//             //std::cout << e.first << ':' << e.second.first << ' ' << e.second.second << std::endl;
        }
        for (int i = 0; i < persistence_pairs_all_dimensions.size(); ++i) {
            //std::cout << "Dimension " << i << std::endl;
            for (const auto &e: persistence_pairs_all_dimensions[i]) {
                //std::cout << i << ": " << e.first << " -> " << e.second << std::endl;
            }
            //std::cout << std::endl;
//             //std::cout << e.first << ':' << e.second.first << ' ' << e.second.second << std::endl;
        }

        //std::cout << "have done this????" << std::endl;


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



//        phat::boundary_matrix<phat::vector_vector> boundary_matrix;
//        boundary_matrix.set_num_cols(total_size);

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
//            if (simplex_map.find(current_bigger_simplex.second) == simplex_map.end()) {
//                exit(1);
//            }
            auto current_bigger_simplex_second = current_bigger_simplex->get_simplex_by_last_node();
//            int row_number = bigger_simplex_map[current_bigger_simplex];
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

        std::ofstream fout("/Users/leonardbee/Desktop/boundary.txt");
//        fout << "start\n";
        fout << from_index_to_birth.size() << '\n';
        for (int i = 0; i < from_index_to_birth.size(); ++i) {
            fout << i << ' ' << from_index_to_birth[i] << ' ' << dimensions[i] << '\n';
        }
//        fout << "-1\n";
        for (int i = 0; i < new_matrix.size(); ++i) {
            fout << i << ' ' << new_matrix[i].size() << ' ';
            for (const auto &e: new_matrix[i]) {
                fout << e << ' ';
            }
            fout << std::endl;
        }
    }

    void initialize_simplex_tree(int n) {

//        for (const auto &point: subcloud_of_points) {
//            root->insert(std::vector<int>({point}), 0);
//        }
    }

    void insert_all_simplices_including(int number_of_vertices) {
        for (int i = 3; i <= number_of_vertices; ++i) {
//             //std::cout << "started inserting" << i << std::endl;
//            root->print();
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
                //what to do with infinite intevals?!
                if (persistence_pair.first == INT_MAX) {
                    Simplex_tree_node *stn = ((root->simplices[i])[persistence_pair.second]);
                    birth_time = stn->get_birth_time_protected();
                    persistence_diagram[persistence_diagram.size() - 1].emplace_back(std::make_pair(birth_time, max_radii));
                    death_time = birth_time;
                } else {
                    if (persistence_pair.second > root->simplices[i - 1].size() - 1) {

//                        number_of_errors++;
                        continue;
                    }
//
//
                    auto found_simplex = (root->simplices[i - 1])[root->simplices[i - 1].size() - 1 -
                                                                  persistence_pair.second];
//                for (const auto& e: found_simplex) {
//
//                }
//
                    auto node_with_this_simplex = (((root->simplices[i - 1])[root->simplices[i - 1].size() - 1 -
                                                                             persistence_pair.second]));
                    birth_time = node_with_this_simplex->get_birth_time_protected();
                }
                if (persistence_pair.first == INT_MAX) {
                } else {
                    if (persistence_pair.first < 0 || persistence_pair.first >= (root->simplices[i]).size()) {

                    } else {
                        death_time = (((root->simplices[i])[persistence_pair.first]))->get_birth_time_protected();
                    }
                }
                if (abs(death_time - birth_time) > epsilon_a) {
                    persistence_intervals.emplace_back(std::make_pair(birth_time, death_time));
//

                } else {

//
                }
//
            }
//
//            persistence_diagram.emplace_back(persistence_intervals);

            if (i == 1) {
//                persistence_intervals = {};
                for (int ip = 0; ip < already_paired_bigger_simplices.size(); ++ip) {
                    if (!already_paired_bigger_simplices[ip]) {
                        persistence_intervals.emplace_back(std::make_pair(
                                ((root->simplices[i - 1])[root->simplices[i - 1].size() - 1 -
                                                          ip])->get_birth_time_protected(),
                                max_radii));
                    }
                }

            }
            persistence_diagram.emplace_back(persistence_intervals);


        }
        return persistence_diagram;
    }

    SubCloud(const Cloud &cloud, const std::vector<int> &points, double max_radius, int max_number_of_points_in_simplex) {
        max_radii = max_radius;
        matrix_of_distances = &cloud.distances;
        subcloud_of_points = points;
        root = new Simplex_tree(cloud.simplex_tree, cloud.simplex_tree->simplices.size());

//        for (int i = 0; i < subcloud_of_points.size(); ++i) {
//
//        }
        for (const auto &point: points) {
            auto found = root->find(std::vector<int>({point}), -1);
            root->simplices[0].emplace_back(found);
        }


        adjacency_matrix = std::vector<boost::dynamic_bitset<>>(points.size(), boost::dynamic_bitset<>(points.size()));
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
//                        if (points[i] == 24 && points[j] == 27)
//                            exit(3);
                        root->insert(std::vector<int>({points[i], points[j]}),
                                     cloud.distances[std::max(points[i], points[j])][std::min(points[i], points[j])]);
                    }
                }
            }
        }
//        root->print();
    }

    friend std::ostream &operator<<(std::ostream &os, const std::vector <boost::dynamic_bitset<>> &matrix);
};


std::ostream &operator<<(std::ostream &os, const SubCloud &subcloud) {
    os << "Printing SubCloud" << std::endl;
    for (int i = 0; i < subcloud.adjacency_matrix.size(); ++i) {
        os << std::setw(space) << subcloud.subcloud_of_points[i] << ": ";
        for (int j = 0; j < subcloud.adjacency_matrix[i].size(); ++j) {
            os << subcloud.adjacency_matrix[i][j];
        }
        os << std::endl;
    }
    return os;
}


std::ostream &operator<<(std::ostream &os, const std::vector <std::vector<double>> &matrix) {
//    os << "Printing distances" << std::endl;
    for (const auto &row: matrix) {
        for (const auto &element: row) {
            os << element << " ";
        }
        os << std::endl;
    }
    return os;
}
//
//template<typename T>
//std::ostream &operator<<(std::ostream &os, const std::vector <T> &matrix) {
//    for (const auto &row: matrix) {
//        os << std::setw(space) << row;
//    }
//
//    return os;
//}





//template<typename T>
//std::ostream &operator<<(std::ostream &os, const std::vector<T> &matrix) {
////    os << "Printing distances" << std::endl;
//    for (const auto &row: matrix) {
//        os << std::setw(space) << row;
//    }
//    os << std::endl;
//    return os;
//}



//std::vector <std::vector<std::pair < double, double>>>
//
//get_persistence_pairs(const Cloud &cloud, int count) {
//    std::vector<int> dots;
//
//    int max_number_of_points_in_simplex = 4;
//    int number_of_dots = number_of_dots_in_file;
//    std::set<int> not_include = {0, 2, 6, 7, 8, 9, 10, 23, 25};
//    std::set<int> nnot_include = {0, 1, 2, 6, 8, 23, 24, 25};
//    if (count == 0) {
//        for (int i = 0; i < number_of_dots; ++i) {
//            dots.emplace_back(i);
//        }
//    } else if (count == 1) {
//        for (int i = 0; i < number_of_dots; ++i) {
//            if (not_include.find(i) == not_include.end()) {
//                dots.emplace_back(i);
//            }
//        }
//    } else {
//        for (int i = 0; i < number_of_dots; ++i) {
//            if (nnot_include.find(i) == nnot_include.end()) {
//                dots.emplace_back(i);
//            }
//        }
//    }
//
//    auto subcloud = new SubCloud(cloud, dots, 8, max_number_of_points_in_simplex);
//    {
//        mute.lock();
//        std::cout << "started building tree" << std::endl;
//        mute.unlock();
//    }
//    subcloud->insert_all_simplices_including(max_number_of_points_in_simplex);
//    std::vector < std::vector < std::pair < double, double >> > result;
//    {
//        mute.lock();
//        std::cout << "started reducing boundary matrix compressed" << std::endl;
//        mute.unlock();
//    }
//    subcloud->get_all_dimensions_landscape<phat::chunk_reduction, phat::vector_vector>(); //58 Mb memory 332 (282) ms (50 dots r = 8)
//
////    subcloud->get_all_dimensions_landscape<phat::twist_reduction, phat::vector_vector>(); //56 Mb memory 277 (222) ms (50 dots r = 8)
//
////    subcloud->get_all_dimensions_landscape<phat::row_reduction, phat::vector_vector>>(); //108 Mb memory 452 ms (50 dots r = 8)
//    delete subcloud->root;
//    delete subcloud;
//    return result;
//}

std::vector<int> get_random_sample(std::vector<int> &vector_of_points, int size_of_one_sample) {
    std::random_device rd;
    std::mt19937 g(rd());
    shuffle(vector_of_points.begin(), vector_of_points.end(), g);
    for (const auto &e: vector_of_points) {
        //std::cout << e << " ";
    }
    //std::cout << std::endl;
    sort(vector_of_points.begin(), vector_of_points.begin() + size_of_one_sample);
    std::vector<int> subcloud(vector_of_points.begin(), vector_of_points.begin() + size_of_one_sample);
    std::vector<int> set_of_points = std::vector<int>(vector_of_points.begin(),
                                                      vector_of_points.begin() + size_of_one_sample);
    //std::cout << std::endl;
    for (const auto &e: set_of_points) {
        //std::cout << e << " ";
    }
    //std::cout << std::endl;
//    m_cntr.lock();
//    auto subset = vector_of_sets[cntr];
//    ++cntr;
//    m_cntr.unlock();
    return set_of_points;
}

void get_persistence_pairs_sparse(Cloud* cloud, double radii, double subsample_density_coefficient,
        tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> &v_pairs) {



    std::cout << "non flag entered pairs construction" << std::endl;

    if (flag) {
        mute.lock();
        std::cout << "entered pairs construction" << std::endl;
        mute.unlock();
    }
    int max_number_of_points_in_simplex = 4;
    int number_of_dots = cloud->size;
    std::vector<int> dots(number_of_dots);

    for (int i = 0; i < number_of_dots; ++i) {
        dots[i] = i;
    }

    auto subsample = get_random_sample(dots, (int) ((dots.size()) * subsample_density_coefficient));
    if (flag) {
        mute.lock();
        std::cout << "trying construct tree" << std::endl;
        mute.unlock();
    }
    auto subcloud = new SubCloud(*cloud, subsample, radii, max_number_of_points_in_simplex);

    if (flag) {
        mute.lock();
        std::cout << "constructed tree" << std::endl;
        mute.unlock();
    }


    std::cout << number_of_dots << " started building tree \n\n\n\n" << subsample.size() <<std::endl;
    for (const auto &e: subsample) {
        std::cout << e << ',';
    }
    std::cout << "_________________________________________________________________________" << std::endl;


    subcloud->insert_all_simplices_including(max_number_of_points_in_simplex);
    std::vector <std::vector<std::pair < double, double>>>result;

    mute.lock();

    int ss = 0;
    for (const auto &e: (subcloud->root->simplices)) {
        ss += e.size();
    }

    std::cout << "started calculating boundary matrix compressed " << ss <<std::endl;
    matrix_size_cntr += ss;
    mute.unlock();
    //only  1



    auto persistence_landscape = subcloud->get_all_dimensions_landscape(); //flag cohomology //58 Mb memory 332 (282) ms (50 dots r = 8)

    std::cout << "delete subcloud root" << std::endl;
//    delete subcloud->root;
    std::cout << "delete subcloud" << std::endl;
    delete subcloud;

    mute.lock();
    for (int i = 0; i < persistence_landscape.size(); ++i) {
        std::cout << "DIM " << i << std::endl;
        sort(persistence_landscape[i].begin(), persistence_landscape[i].end());
        for (const auto& e: persistence_landscape[i]) {
            std::cout << e.first << " " << e.second << std::endl;
        }
    }
    v_pairs.push_back(persistence_landscape);
    mute.unlock();

}

void
get_average_landscape_once(tbb::concurrent_vector<std::vector<std::pair<double, double>>>& diagram,
        Cloud* cloud, int number_of_thread_workers, double radii = 0.5,
                           double subsample_density_coefficient = 0.3,
                           int number_of_samples = 10, bool print_pairs = false) {

    boost::asio::thread_pool pool(number_of_thread_workers);
    if (flag) {
        std::cout << "Created pool" << std::endl;
    }
    tbb::concurrent_vector < tbb::concurrent_vector < std::vector < std::pair < double, double >> >>
                                                                                               all_persistence_diagrams;

    for (int i = 0; i < number_of_samples; ++i) {
        boost::asio::post(pool,
                          bind(get_persistence_pairs_sparse, cloud, radii, subsample_density_coefficient,
                                    ref(all_persistence_diagrams)));
    }
    pool.join();
    {
        mute.lock();
        std::cout << "joined" << std::endl;
        mute.unlock();
    }
    int counter = 1;
    std::vector <std::vector<Persistence_landscape>> persistence_landscapes(cloud->dimension);
    if (all_persistence_diagrams.empty()) {
        return;
    }
    diagram = all_persistence_diagrams[0];
    get_average_landscape(all_persistence_diagrams);
}

double main_algorithm(tbb::concurrent_vector<std::vector<std::pair<double, double>>>& diagram,
                        std::string from, std::string to,
                        int max_rank,
                        double max_edge_length,
                        bool gudhi_format,
                        int number_of_thread_workers = 1,
                        int number_of_samples = 1,
                        double subsample_density_coefficient = 1.0,
                        bool print_pairs = false) {
    zero_cntr = 0;
    extra_cntr = 0;
    matrix_size_cntr = 0;
    flag = true;
    number_of_dots_in_file = 0;
//    cntr = 0;
    int max_number_of_points_in_simplex = max_rank + 2;
    if (flag) {
        std::cout << "Trying to create" << std::endl;
    }
    Cloud *matrix = new Cloud(from, max_number_of_points_in_simplex);
    if (flag) {
        std::cout << "Created" << std::endl;
    }

//    std::ofstream fout("output.txt");
//    fout << matrix->distances;
    std::vector<int> dots;

    int number_of_dots = matrix->size;
    std::set<int> not_include = {0, 2, 6, 8, 23, 25};
    for (int i = 0; i < number_of_dots; ++i) {
        if (not_include.find(i) == not_include.end()) {
            dots.emplace_back(i);
        }
    }

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    if (flag) {
        std::cout << "Started calculating" << std::endl;
    }
    get_average_landscape_once(diagram, matrix, number_of_thread_workers, max_edge_length, subsample_density_coefficient,
                          number_of_samples);
    if (flag) {
        std::cout << "Finished calculating" << std::endl;
    }
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();


    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
//    delete matrix->simplex_tree;
    std::cout << "Deleting tree" << std::endl;
    delete matrix;

    std::cout << zero_cntr << " so " << extra_cntr << std::endl << "duration " << duration << std::endl;

    return duration;
}


