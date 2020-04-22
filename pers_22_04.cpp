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
#include <boost/dynamic_bitset.hpp>
#include <tbb/concurrent_unordered_map.h>
#include <boost/container/flat_map.hpp>
#include <unordered_map>

#define space 5
#define max_rad 100

#define epsilon (1e-10)

using namespace std;


mutex mute;

int number_of_dots_in_file;
int number_of_errors = 0;

ostream &operator<<(ostream &os, const vector<boost::dynamic_bitset<>> &matrix) {
    os << endl;
    for (int j = 0; j < matrix.size(); ++j) {
        boost::dynamic_bitset<> row = matrix[j];
        os << setw(space) << j << ": ";
        for (int i = 0; i < row.size(); ++i) {
            os << row[i];
        }
        os << endl;
    }
    return os;
}

void _xor_vector(vector<int> &l, const vector<int> &r) {
    vector<int> resulting_vector;
    resulting_vector.reserve(l.size());
    int l_pointer = 0, r_pointer = 0;
    while (true) {
        if (l_pointer >= l.size() && r_pointer >= r.size()) {
            break;
        }
        int l_temp = INT_MAX, r_temp = INT_MAX;
        if (l_pointer < l.size()) {
            l_temp = l[l_pointer];
        }
        if (r_pointer < r.size()) {
            r_temp = r[r_pointer];
        }
        bool left_is_less = (l_temp < r_temp);
        bool are_equal = (l_temp == r_temp);
        if (are_equal) {
            ++l_pointer;
            ++r_pointer;
            continue;
        }
        if (left_is_less) {
            resulting_vector.emplace_back(l_temp);
            ++l_pointer;
        } else {
            resulting_vector.emplace_back(r_temp);
            ++r_pointer;
        }
    }
    swap(l,resulting_vector);
}
void xor_vector(vector<int> &l, vector<int> &r) {
    std::vector<int> v(l.size() + r.size());                      // 0  0  0  0  0  0  0  0  0  0
    std::vector<int>::iterator it;

    it = std::set_symmetric_difference(l.begin(), l.end(), r.begin(), r.end(), v.begin());
    v.resize(it - v.begin());
//    std::vector<int>
    auto v1(v);
//    l = std::move(v);
//    auto vector_copy(std::move(v));
    std::swap(l, v1);
}


struct Simplex_tree_node {
protected:
    int vertex_number;
    float birth_time = - epsilon;
    Simplex_tree_node *previous_vertex = nullptr;

public:
    float get_signed_birth_time() {
        return birth_time;
    }
    int get_vertex_number() {
        return vertex_number;
    }
    void set_previous(Simplex_tree_node *simplex_tree_node) {
        previous_vertex = simplex_tree_node;
    }

    float get_birth_time_protected() {
        return (abs(birth_time) < 2 * epsilon ) ? (0) : abs(birth_time);
    }

    float get_birth_time() {
        return (abs(birth_time) < 2 * epsilon ) ? (0) : abs(birth_time);
    }


    void set_prev(Simplex_tree_node *prev_node) {
        previous_vertex = prev_node;
    }


//    vector<int> get_simplex_by_last_node() {
//        auto current_node = this;
//        vector<int> result;
//        while (current_node != nullptr) {
//            result.emplace_back(current_node->vertex_number);
////            //cout << "\n->" << current_node->vertex_number << endl;
//            current_node->previous_vertex->sh_mutex.lock_shared();
//            auto prev = current_node;
//            current_node = current_node->previous_vertex;
//            prev->sh_mutex.unlock_shared();
//        }
//        reverse(result.begin(), result.end());
//        return result;
//    }

    vector<int> get_simplex_by_last_node(int size) {
        auto current_node = this;
        vector<int> result(size);
        int p = size - 1;
        while (current_node != nullptr) {
            result[p] = abs(current_node->vertex_number);
            --p;
            current_node = current_node->previous_vertex;

        }

        return result;
    }

    Simplex_tree_node() {

    }

    Simplex_tree_node(int current_vertex_number, float current_birth_time) :
            vertex_number(current_vertex_number),
            birth_time(-current_birth_time) {
//        cout << "passed " << current_vertex_number << " vs " << vertex_number << endl;
    }

};




struct Simplex_tree_node_inner: Simplex_tree_node {
private:
    mutable shared_mutex sh_mutex;
//    mutable mutex sh_mutex;
//    std::unordered_map<int, Simplex_tree_node *> successors;
//    unordered_map<int, Simplex_tree_node *> successors;
//    tbb::concurrent_unordered_map<int, Simplex_tree_node *> successors;
//    unordered_map<int, Simplex_tree_node *> successors;
    boost::container::flat_map<int, Simplex_tree_node *> successors;
public:
//    void set_previous(Simplex_tree_node *simplex_tree_node) {
//        previous_vertex = simplex_tree_node;
//    }
//
//    float get_birth_time_protected() {
//        float res = birth_time;
//        return res;
//    }
//
//    float get_birth_time() {
//        return birth_time;
//    }

    void lock_unique() {
        sh_mutex.lock();
    }

    void unlock_unique() {
        sh_mutex.unlock();
    }

    void lock_shared() {
        sh_mutex.lock_shared();
//        sh_mutex.lock();
    }

    void unlock_shared() {
        sh_mutex.unlock_shared();
//        sh_mutex.unlock();
    }

//    void set_prev(Simplex_tree_node *prev_node) {
//        previous_vertex = prev_node;
//    }

    auto *get_successors() {
        return &successors;
    }

//    vector<int> get_simplex_by_last_node() {
//        auto current_node = this;
//        vector<int> result;
//        while (current_node != nullptr) {
//            result.emplace_back(current_node->vertex_number);
////            //cout << "\n->" << current_node->vertex_number << endl;
//            current_node->previous_vertex->sh_mutex.lock_shared();
//            auto prev = current_node;
//            current_node = current_node->previous_vertex;
//            prev->sh_mutex.unlock_shared();
//        }
//        reverse(result.begin(), result.end());
//        return result;
//    }

//    vector<int> get_simplex_by_last_node(int size) {
//        auto current_node = this;
//        vector<int> result(size);
//        int p = size - 1;
////        lock_shared();
//        while (current_node != nullptr) {
//            result[p] = current_node->vertex_number;
//            --p;
//            //changed
////            current_node->unlock_shared();
////            if (current_node->previous_vertex) {
////                current_node->previous_vertex->lock_shared();
////            }
//            current_node = current_node->previous_vertex;
//
//        }
//
//        return result;
//    }

//    Simplex_tree_node_end() {
//
//    }

//    Simplex_tree_node_end(int current_vertex_number, float current_birth_time) :
//            vertex_number(current_vertex_number),
//            birth_time(current_birth_time) {}

    Simplex_tree_node_inner() {
    }
    Simplex_tree_node_inner(int current_vertex_number, float current_birth_time) {
//        cout << "constructor called on " << current_vertex_number << " & " << current_birth_time << endl;
//        Simplex_tree_node(current_vertex_number, -current_birth_time);
        vertex_number = current_vertex_number;
        birth_time = current_birth_time;

//        cout << "but got  " << vertex_number << " & " << birth_time << endl << endl;
    }

    ~Simplex_tree_node_inner() {
//        cout << "deleting " << vertex_number << " " << birth_time << endl;
        for (const auto &e: successors) {
            if (e.second) {
                if (e.second->get_signed_birth_time() < 0 && abs(e.second->get_birth_time()) > 2 * epsilon) {
//                    cout << "deleting as final node " << e.second->get_vertex_number() << " " << e.second->get_birth_time() << endl;
                    delete (e.second);
                } else {
//                    cout << "deleting as inner node " << e.second->get_vertex_number() << " " << e.second->get_birth_time() << endl;
                    delete ((Simplex_tree_node_inner*)e.second);
                }
            }
        }
    }

    void print(int level) {
        for (int i = 0; i < level; ++i) {
            cout << "*";
        }
        cout << " " << get_vertex_number() << " [" << this->get_birth_time() << "]" << endl;
        for (const auto &e: successors) {
            if (e.second->get_signed_birth_time() < 0) {
                for (int i = 0; i <= level; ++i) {
                    cout << "*";
                }
                cout << " " << e.second->get_vertex_number() << " [" << e.second->get_signed_birth_time() << "]" << endl;
            } else {
                ((Simplex_tree_node_inner *) e.second)->print(level + 1);
            }
        }
    }

    void print() {
        cout << endl;
        print(1);
    }
};

struct Simplex_tree {     //this is only 1's thread property
private:
    Simplex_tree_node_inner *all_first_vertices; //shared resource
public:
    bool to_delete = true;
    vector<vector<Simplex_tree_node *>> simplices; //private for thread-owner (now not but should be!!!)

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
    Simplex_tree(int n, int cloud_size) {
        all_first_vertices = new Simplex_tree_node_inner();
        simplices = vector<vector<Simplex_tree_node *>>(n);
        for (int i = 0; i < cloud_size; ++i) {
            auto new_node = new Simplex_tree_node_inner(i, epsilon);
            (*(((Simplex_tree_node_inner*)all_first_vertices)->get_successors()))[i] = new_node;
            simplices[0].emplace_back(new_node);
        }

//        all_first_vertices->print();
//        cout << "constructor done " << endl;
    }

    Simplex_tree(Simplex_tree *simplex_tree, int n) {
        to_delete = false;
        all_first_vertices = simplex_tree->all_first_vertices;
        simplices = vector<vector<Simplex_tree_node *>>(n);
    }


    Simplex_tree_node *find(const vector<int> &simplex, int not_included) {
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
            if (i != 0 && i < simplex.size() - 1) {
                ((Simplex_tree_node_inner*)current_simplex_tree_node)->lock_shared();
            }

            Simplex_tree_node* next_node = (*(((Simplex_tree_node_inner*)current_simplex_tree_node)->get_successors()))[simplex[i]];

            if (i != 0 && i < simplex.size() - 1) {
                ((Simplex_tree_node_inner*)current_simplex_tree_node)->unlock_shared();
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

    void insert(Simplex_tree_node *current_simplex_tree_node_temp, int i, float max_birth_time, int vertex) {
        bool is_last = (max_birth_time < 0);
        max_birth_time = abs(max_birth_time);

        auto current_simplex_tree_node = (Simplex_tree_node_inner*)current_simplex_tree_node_temp;
        if (!is_last) {

        }
        current_simplex_tree_node->lock_unique();
        auto found = current_simplex_tree_node->get_successors()->find(vertex);

        auto new_simplex_node = new Simplex_tree_node(vertex, max_birth_time);
//        cout << " here looking for " << endl;
//        if (current_simplex_tree_node->get_successors()->find(vertex) ==
//            current_simplex_tree_node->get_successors()->end()) {
            if ((current_simplex_tree_node->get_successors()->insert(make_pair(vertex, new_simplex_node))).second) {

//            cout << "not found " << endl;

            if (is_last) {

                if (i > 0) {
                    new_simplex_node->set_previous(current_simplex_tree_node);
                }
                (*current_simplex_tree_node->get_successors())[vertex] = new_simplex_node;

                simplices[i].emplace_back(new_simplex_node);
            } else {
                auto new_simplex_node = new Simplex_tree_node_inner(vertex, max_birth_time);

                if (i > 0) {
                    new_simplex_node->set_previous(current_simplex_tree_node);
                }
                (*current_simplex_tree_node->get_successors())[vertex] = new_simplex_node;

                simplices[i].emplace_back(new_simplex_node);
            }
            //}
        } else {
//            cout << "found" << std::endl;
            auto found = current_simplex_tree_node->get_successors()->find(vertex);
            simplices[i].emplace_back(found->second);
            delete new_simplex_node;
        }


//        auto next = (*current_simplex_tree_node->get_successors())[vertex];
//        if (next && i > 0) {
//            next->previous_vertex = current_simplex_tree_node;
//        }
        current_simplex_tree_node->unlock_unique();
//        next->unlock_unique();
    }

//    Simplex_tree_node *find(const vector<int> &simplex) {
//        Simplex_tree_node *current_simplex_tree_node = all_first_vertices;
//        if (current_simplex_tree_node == nullptr) {
//            return nullptr;
//        }
//        for (int i = 0; i < simplex.size(); ++i) {
//            current_simplex_tree_node = current_simplex_tree_node->successors[simplex[i]];
//
//            if (current_simplex_tree_node == nullptr) {
//                return nullptr;
//            }
//        }
//        return current_simplex_tree_node;
//    }

    void insert(const vector<int> &simplex, float birth_time) {

        auto node = find(simplex, simplex.size() - 1);
        insert(node, simplex.size() - 1, birth_time, simplex[simplex.size() - 1]);
//        auto size = simplex.size();
//        Simplex_tree_node *current_simplex_tree_node = all_first_vertices;
//        Simplex_tree_node *previous_simplex_tree_node = nullptr;
//
//        if (!current_simplex_tree_node) {
//            return;
//        }
//        if (size == 0) {
//            return;
//        }
//        current_simplex_tree_node->lock_unique();
//        for (int i = 0; i < simplex.size(); ++i) {
//            float max_birth_time = birth_time;
//            if (i == simplex.size() - 1) {
//                auto new_simplex_node = new Simplex_tree_node(simplex[i], max_birth_time);
//                if (i == simplex.size() - 1) {
//                    current_simplex_tree_node->get_successors()->insert({simplex[i], new_simplex_node});
//                }
//                simplices[i].emplace_back(make_pair(max_birth_time, new_simplex_node));
//            }
//
//            auto next = (*current_simplex_tree_node->get_successors())[simplex[i]];
//            next->lock_unique();
//
//            if (next && i > 0) {
//                next->previous_vertex = current_simplex_tree_node;
//            }
//            current_simplex_tree_node->unlock_unique();
//
//            current_simplex_tree_node = next;
//        }
//        current_simplex_tree_node->unlock_unique();

    }

    void insert_all_simplices(const vector<boost::dynamic_bitset<>> &matrix_of_adjacency,
                              unordered_map<int, int> &new_order_of_points,
                              const vector<int> &subcloud_of_points,
                              int number_of_vertices,
                              const vector<vector<float>> *matrix_of_distances) {
        if (number_of_vertices < 2 || matrix_of_adjacency.empty()) {
            return;
        }
        if (number_of_vertices - 2 >= simplices.size()) {
            return;
        }
        //cout << "subcloud: ";
//        for (const auto &e: subcloud_of_points) {
//            //cout << e << " ";
//        }
        //cout << endl;
        //cout << "less simplices-vector size is " << simplices[number_of_vertices - 2].size() << endl;
        boost::dynamic_bitset<> row(matrix_of_adjacency[0].size());
        //cout << matrix_of_adjacency;
        for (const auto &temp_subsimplex: simplices[number_of_vertices - 2]) {


            row.set();
            int max_vertex = INT_MIN;

            auto temp_copy_of_subsimplex = temp_subsimplex->get_simplex_by_last_node(number_of_vertices - 1);
//            for (const auto& e: temp_copy_of_subsimplex) {
//                //cout << e << " ";
//            }
//            //cout << endl;
//
//            if (temp_subsimplex2.second.size() > number_of_vertices - 1)
//                exit(17);
            for (const auto &vertex: temp_copy_of_subsimplex) {
                row &= matrix_of_adjacency[new_order_of_points[vertex]];
                max_vertex = max(new_order_of_points[vertex], max_vertex);
            }
            if (max_vertex != INT_MIN) {
                auto temp = row.find_next(max_vertex);
                while (temp < row.size()) {
//                    temp_copy_of_subsimplex = temp_subsimplex2.second;
                    auto found_node = find(temp_copy_of_subsimplex, -1);
                    float max_birth_time = found_node->get_birth_time_protected();

                    for (int i = 0; i < temp_copy_of_subsimplex.size(); ++i) {
                        max_birth_time = max(max_birth_time,
                                             (*matrix_of_distances)[max(temp_copy_of_subsimplex[i],
                                                                        subcloud_of_points[temp])]
                                             [min(temp_copy_of_subsimplex[i], subcloud_of_points[temp])]);
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
        //cout << "inserted bigger simplices: " << endl;
        //cout << simplices[number_of_vertices - 1].size() << endl;
    }

//    void print() {
//        (all_first_vertices)->print();
//    }
    void print() {
        for (const auto &e: *(((Simplex_tree_node_inner*)all_first_vertices)->get_successors())) {
            //cout << endl << e.first << ": " << endl;
            if (e.second) {
                ((Simplex_tree_node_inner*)e.second)->print();
            }
        }
//        for (const auto &simplices_array: simplices) {
//            for (const auto &simplex: simplices_array) {
//                //cout << simplex.first << ": ";
//                for (const auto &e: simplex.second->get_simplex_by_last_node()) {
//                    //cout << setw(space) << e;
//                }
//                //cout << endl;
//            }
//        }
    }

    ~Simplex_tree() {
        if (to_delete) {
//            cout << "ready to delete root" << endl;
            delete ((Simplex_tree_node_inner*)all_first_vertices);
//            cout << "deleted root" << endl;

            all_first_vertices = nullptr;
        }
    }

};


struct Cloud {
    unsigned int dimension, size;
    vector<vector<float>> points;
    vector<vector<float>> distances;
    Simplex_tree *simplex_tree;

    Cloud(const string &path, int n) {
        std::ifstream fin(path);
        if (!fin.is_open()) {
            cout << "problems opening file" << endl;
            exit(1);
        }
        fin >> size >> dimension;
        simplex_tree = new Simplex_tree(n, size);
//        simplex_tree->print();
//        cout << "printed|n\n\n\n" << endl;
        number_of_dots_in_file = size;

        points = vector<vector<float>>(size, vector<float>(dimension, 0));
        for (int i = 0; i < points.size(); ++i) {
            for (int j = 0; j < dimension; ++j) {
                fin >> points[i][j];
            }
        }
        distances = vector<vector<float>>(size);
        for (int i = 0; i < distances.size(); ++i) {
            distances[i] = vector<float>(i);
            for (int j = 0; j < i; ++j) {
                float sum_of_squared_differences = 0;
                for (int k = 0; k < dimension; ++k) {
                    sum_of_squared_differences += (points[i][k] - points[j][k]) * (points[i][k] - points[j][k]);
                }
                distances[i][j] = sqrt(sum_of_squared_differences);
            }
        }
//        for (int i = 0; i < points.size(); ++i) {
//            simplex_tree->insert(vector<int>({i}), 0);
//        }
//        for (const auto &row: distances) {
//            for (const auto &e: row) {
//                cout << setw(space + 5) << e;
//            }
//            cout << endl;
//        }

    }

    ~Cloud() {

    }

    friend ostream &operator<<(ostream &os, const vector<vector<float>> &matrix);

    void print_distances(ostream &os) {
        //cout << "Printing distances" << endl;
        for (const auto &row: distances) {
            for (const auto &element: row) {
                //cout << setw(space) << element;
            }
            //cout << endl;
        }
    }
};




struct Betti_matrix {
    int number_of_points_in_the_less_simplex;
    Simplex_tree *simplex_tree;
    vector<Simplex_tree_node *> *less_simplices;
    vector<Simplex_tree_node *> *bigger_simplices;
    vector<int> *original_numbers_of_points;

    unordered_map<Simplex_tree_node *, int>& less_simplex_map_;
//    unordered_map<Simplex_tree_node *, int>& bigger_simplex_map;

    Betti_matrix(int new_number_of_points_in_the_less_simplex,
                 vector<Simplex_tree_node *> *new_less_simplices,
                 vector<Simplex_tree_node *> *new_bigger_simplices,
                 Simplex_tree *new_simplex_tree,
                 vector<int> *new_original_numbers_of_points,
                 unordered_map<Simplex_tree_node *, int>& new_less_simplex_map,
                 unordered_map<Simplex_tree_node *, int>& new_bigger_simplex_map) :
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
    vector<vector<int>>* new_matrix_pointer = nullptr;



    vector<pair<int, int>> reduce_matrix(vector<vector<int>> &temp_betti_matrix,
                                         boost::dynamic_bitset<> &paired_bigger_simplices_should_not_be_used,
                                         bool is_max_simplex_rank,
                                         priority_queue<int, vector<int>, function<bool(const int&, const int&)>> &rows_with_lowest_bits) {
        int xor_counter = 0;
        int max_size = 0;
        int mm_size = 0;
        for (const auto& e: temp_betti_matrix) {
            max_size += e.size();
        }

        //cout << "rows " << rows_with_lowest_bits.size() << endl;
        vector<pair<int, int>> persistence_pairs;
        int max_number_of_vertices1 = INT_MIN, max_number_of_vertices2 = INT_MIN;
        persistence_pairs.reserve(min(less_simplices->size(), bigger_simplices->size()));
        int pos_in_paired_simplices = paired_bigger_simplices_should_not_be_used.find_first();    //this is reversed vector
//        int pos_in_paired_simplices = paired_bigger_simplices_should_be_used.find_first();
        boost::dynamic_bitset<> should_use(temp_betti_matrix.size());
//        boost::dynamic_bitset<> paired_simplices(betti_matrix.size());
//        paired_simplices.set();
        should_use.set(); //while normalizing matrix we have to mark used rows in order not to normalize them again
        if (!is_max_simplex_rank) {
            while (pos_in_paired_simplices < paired_bigger_simplices_should_not_be_used.size()) {
                ////            betti_matrix[pos_in_paired_simplices].reset();
                temp_betti_matrix[pos_in_paired_simplices] = {};
                should_use[pos_in_paired_simplices] = 0;
                pos_in_paired_simplices = paired_bigger_simplices_should_not_be_used.find_next(pos_in_paired_simplices);
            }
            //cout << "isnt max 1" << endl;
        }
        boost::dynamic_bitset<> temp_paired_simplices;
        if (!is_max_simplex_rank) {
            temp_paired_simplices = boost::dynamic_bitset<>(paired_bigger_simplices_should_not_be_used);
        }
        paired_bigger_simplices_should_not_be_used = boost::dynamic_bitset<>(less_simplices->size(), 0); //we will pass this bitset consisting of already paired simplices to next iteration
        //but it has to be resized (to size of less simplices)
//        paired_bigger_simplices_should_be_used.set();

        if (temp_betti_matrix.size() == 0) {
            return persistence_pairs;
        }
//        pair<int, int> lowest_bit_and_row = {INT_MAX, INT_MAX};
        int lowest_row = INT_MAX;
        //cout << "got here" << endl;

        int max_counter = INT_MIN;
        while (!rows_with_lowest_bits.empty()) {
            //cout << "new size " << rows_with_lowest_bits.size() << endl;
//            lowest_bit_and_row = rows_with_lowest_bits.top();
            lowest_row = rows_with_lowest_bits.top();
            int lowest_bit = temp_betti_matrix[lowest_row][0];
            rows_with_lowest_bits.pop();
            if (should_use[lowest_row]) {
                should_use[lowest_row] = 0;
                persistence_pairs.emplace_back(make_pair(lowest_row, lowest_bit));
                paired_bigger_simplices_should_not_be_used[paired_bigger_simplices_should_not_be_used.size() - 1 -
                        lowest_bit] = 1;
            }
            if (rows_with_lowest_bits.empty()) {
                break;
            }
            auto next_top_bit_and_row = rows_with_lowest_bits.top();
            while (!rows_with_lowest_bits.empty() && temp_betti_matrix[rows_with_lowest_bits.top()][0]/*rows_with_lowest_bits.top().first*/ == lowest_bit) {
                next_top_bit_and_row = rows_with_lowest_bits.top();
                rows_with_lowest_bits.pop();

                max_size -= temp_betti_matrix[next_top_bit_and_row].size();
                max_size -= temp_betti_matrix[lowest_row].size();
                auto temp_max_counter = std::max(temp_betti_matrix[next_top_bit_and_row].size(), temp_betti_matrix[lowest_row].size());

                xor_vector(temp_betti_matrix[next_top_bit_and_row], temp_betti_matrix[lowest_row]);
                auto temp2_max_counter = std::max(temp_betti_matrix[next_top_bit_and_row].size(), temp_betti_matrix[next_top_bit_and_row].size());

                max_size += temp_betti_matrix[next_top_bit_and_row].size();
                max_size += temp_betti_matrix[lowest_row].size();
                mm_size = std::max(mm_size, max_size);
                max_counter = std::max(max_counter, (int)std::max(temp2_max_counter, temp_max_counter));
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

        }
        //cout << "here" << endl;
        if (!is_max_simplex_rank) {
            for (int i = 0; i < temp_betti_matrix.size(); ++i) {
                if (!is_max_simplex_rank && temp_betti_matrix[i].empty() && !temp_paired_simplices[i]) {
                    persistence_pairs.emplace_back(make_pair(INT_MAX, i));
                }
            }
            //cout << "isnt max 2" << endl;
        }
        mute.lock();
        cout << "xor number is " << xor_counter << endl << "max number in a row is " << max_counter << endl;
        cout << "size is " << mm_size << endl;
        mute.unlock();
        return persistence_pairs;
    }

    vector<pair<int, int>> construct_betti_matrix(boost::dynamic_bitset<> &already_paired_simplices,
                                                  bool is_max_simplex_rank,
                                                  int current_number_of_vertices) {
        bool f = true;
        mute.lock();
        cout << "constructing matrices!" << std::endl;

        mute.unlock();
        int simplex_indicator;

        simplex_indicator = 0;

//        boost::container::flat_map<Simplex_tree_node*, int> less_simplex_map;
        unordered_map<Simplex_tree_node*, int> less_simplex_map;
        less_simplex_map.reserve(less_simplices->size());
        for (const auto &current_less_simplex: *less_simplices) {
//            simplex_map.emplace(make_pair(current_less_simplex, simplex_indicator));

            less_simplex_map.emplace(make_pair(current_less_simplex, simplex_indicator));
            ++simplex_indicator;
        }
//        if (is_max_simplex_rank) {
//            simplex_indicator = 0;
//            bigger_simplex_map.reserve(bigger_simplices->size());
//            for (const auto &current_bigger_simplex: *bigger_simplices) {
//                bigger_simplex_map.emplace(make_pair(current_bigger_simplex, simplex_indicator));
//                ++simplex_indicator;
//            }
//        }


        vector<vector<int>> new_matrix(bigger_simplices->size());

        //cout << "started inserting into priority queue" << endl;
        vector<int> vector_for_pq;


        //cout << " \n\n\n____________________i s max n\n\n\n\n" << endl;
        if (is_max_simplex_rank) {
            vector_for_pq.reserve(bigger_simplices->size());
        }
        new_matrix_pointer = &new_matrix;
        mute.lock();
        cout << "iteration part" << endl;

        cout << "attention" << endl;
        mute.unlock();
        for (int i = 0; i < bigger_simplices->size(); ++i) {
//            if (simplex_map.find(current_bigger_simplex.second) == simplex_map.end()) {
//                exit(1);
//            }
            auto &current_bigger_simplex = (*bigger_simplices)[i];
            auto current_bigger_simplex_second = current_bigger_simplex->get_simplex_by_last_node(
                    current_number_of_vertices + 1);
//            int row_number = bigger_simplex_map[current_bigger_simplex];
            int row_number = i;
//            mute.lock();
//
//            cout << "row number attention!!! " << bigger_simplex_map[current_bigger_simplex] << endl;
//            if ( bigger_simplex_map[current_bigger_simplex] < 0 )
//                exit(51);
//            mute.unlock();
            new_matrix[row_number].reserve(current_bigger_simplex_second.size());
            for (int j = 0; j < current_bigger_simplex_second.size(); ++j) {
                new_matrix[row_number].emplace_back(less_simplices->size() - 1 -
                                                    less_simplex_map[simplex_tree->find(current_bigger_simplex_second, j)]);
            }
            sort(new_matrix[row_number].begin(), new_matrix[row_number].end());
            vector_for_pq.emplace_back(row_number);
        }
        //cout << "finLYY inserting" << endl;

//        priority_queue<pair<int, int>, stack<pair<int, int>>, mycomparison> pq;
        mute.lock();
        cout << "starting priority que size " << 0 << endl;
        cout << "with  " << bigger_simplices->size() << " & " << less_simplices->size() << endl;
        mute.unlock();
//        exit(1);

        priority_queue<int, vector<int>, function<bool(const int&, const int&)>> rows_with_lowest_bits{
                [&](const int& rhs, const int& lhs) -> bool {
//                    return lhs.first > rhs.first;
//                        if (lhs.first > rhs.first) {
//                            return true;
//                        } else {
//                            return lhs.second >= rhs.second;
//                        }
//                    return rhs > lhs;
///////////////smth is wrong - infinite loop
                    return (*new_matrix_pointer)[lhs][0]<(*new_matrix_pointer)[rhs][0] || (!((*new_matrix_pointer)[rhs][0]<(*new_matrix_pointer)[lhs][0]) && lhs<rhs);
//                        return rhs < lhs;
                }, std::move(vector_for_pq)};
        mute.lock();
        cout << "priority que size " << rows_with_lowest_bits.size() << endl;
        mute.unlock();
//        priority_queue<pair<int,int>, vector<pair<int,int>>, function<bool(pair<int,int>, pair<int,int>)>> rows_with_lowest_bits(
//                [&](const pair<int,int>& lhs, const pair<int,int>& rhs) -> bool {
//            return (*new_matrix_pointer)[lhs.second] > (*new_matrix_pointer)[rhs.second];
//            }, vector_for_pq);
//        priority_queue<pair<int, int>, vector<pair<int, int>>, mycomparison> rows_with_lowest_bits (vector_for_pq.begin(), vector_for_pq.end());
//        rows_with_lowest_bits = priority_queue<pair<int, int>, vector<pair<int, int>>, > ;
//        cout << "swapped 0 " << less_simplex_map.size() << " & " << bigger_simplex_map.size() << endl;
//        unordered_map<Simplex_tree_node*, int> new_map_for_simplex_tree_nodes(std::move(bigger_simplex_map));
//        unordered_map<Simplex_tree_node*, int> new_map;
//        std::swap(new_map_for_simplex_tree_nodes, bigger_simplex_map);

//        cout << "swapped 1 " << less_simplex_map.size() << " & " << bigger_simplex_map.size() << endl;
//        std::swap(less_simplex_map, new_map);
//        cout << "swapped 2 " << less_simplex_map.size() << " & " << bigger_simplex_map.size() << endl;
        mute.lock();
        cout << "started reducing" << std::endl;
        mute.unlock();
        return reduce_matrix(new_matrix, already_paired_simplices, is_max_simplex_rank, rows_with_lowest_bits);
    }

};

struct SubCloud {
    unordered_map<int, int> new_order_of_points;
    vector<boost::dynamic_bitset<>> adjacency_matrix;
    vector<int> subcloud_of_points;
    Simplex_tree *root;
    const vector<vector<float>> *matrix_of_distances;
    float max_radii = 1e6;


    unordered_map<Simplex_tree_node *, int> less_simplex_map;
    unordered_map<Simplex_tree_node *, int> bigger_simplex_map;

    ~SubCloud() {

    }

    void initialize_simplex_tree(int n) {

//        for (const auto &point: subcloud_of_points) {
//            root->insert(vector<int>({point}), 0);
//        }
    }

    void insert_all_simplices_including(int number_of_vertices) {
        for (int i = 3; i <= number_of_vertices; ++i) {
//            cout << "started inserting" << i << endl;
//            root->print();
            root->insert_all_simplices(adjacency_matrix, new_order_of_points, subcloud_of_points, i,
                                       matrix_of_distances);
            //cout << "cancelled inserting " << i << endl;
        }
    }

    vector<vector<pair<float, float>>> calculate_betti_matrix(int max_number_of_vertices_in_simplex) {
        bool f = true;
        vector<vector<pair<float, float>>> persistence_diagram;


        boost::dynamic_bitset<> already_paired_bigger_simplices(
                root->simplices[max_number_of_vertices_in_simplex].size(), 0);
        for (int i = 0; i < root->simplices.size(); ++i) {
            auto bigger_simplices = &root->simplices[i];
            stable_sort(bigger_simplices->begin(), bigger_simplices->end(), [](const auto &lhs, const auto &rhs) {
                return lhs->get_birth_time() < rhs->get_birth_time();
            });
        }

        for (int i = max_number_of_vertices_in_simplex; i >= 1; --i) {
            //cout << i + 1 << "-angle into angle-" << i << endl;

            auto betti_matrix = Betti_matrix(i, &root->simplices[i - 1], &root->simplices[i], root,
                                             &subcloud_of_points, less_simplex_map, bigger_simplex_map);


            auto persistence_pairs = betti_matrix.construct_betti_matrix(already_paired_bigger_simplices,
                                                                         i == max_number_of_vertices_in_simplex, i);
            vector<pair<float, float>> persistence_intervals;
            for (const auto &persistence_pair: persistence_pairs) {
                float birth_time = 0;


                float death_time;
                //what to do with infinite intevals?!
                if (persistence_pair.first == INT_MAX) {
                    Simplex_tree_node *stn = ((root->simplices[i])[persistence_pair.second]);
                    birth_time = stn->get_birth_time_protected();
                    persistence_diagram[persistence_diagram.size() - 1].emplace_back(make_pair(birth_time, max_radii));
                    death_time = birth_time;
                } else {
                    if (persistence_pair.second > root->simplices[i - 1].size() - 1) {
                        //cout << "error " << persistence_pair.first << " and " << persistence_pair.second << endl;
                        number_of_errors++;
                        continue;
                    }
//                //cout << "position is " << less_simplices_array.size() - 1 - persistence_pair.second << endl;
//                //cout << less_simplices_array.size() << endl;
                    auto found_simplex = (root->simplices[i - 1])[root->simplices[i - 1].size() - 1 -
                                                                  persistence_pair.second];
//                for (const auto& e: found_simplex) {
//                    //cout << e << " ";
//                }
//                //cout << "\nstarted" << endl;
                    auto node_with_this_simplex = (((root->simplices[i - 1])[root->simplices[i - 1].size() - 1 -
                                                                             persistence_pair.second]));
                    birth_time = node_with_this_simplex->get_birth_time_protected();
                }
                if (persistence_pair.first == INT_MAX) {
                } else {
                    if (persistence_pair.first < 0 || persistence_pair.first >= (root->simplices[i]).size()) {
                        //cout << "ERROR" << endl;
                    } else {
                        death_time = (((root->simplices[i])[persistence_pair.first]))->get_birth_time_protected();
                    }
                }
                if (abs(death_time - birth_time) > epsilon) {
                    persistence_intervals.emplace_back(make_pair(birth_time, death_time));
//                    //cout << persistence_pair.first << " &_ " << persistence_pair.second << endl;

                } else {

//                    //cout << persistence_pair.first << " ±&_ " << persistence_pair.second << endl;
                }
//                //cout << "ended" << endl;
            }
//            //cout << "pushed back" << endl;
//            persistence_diagram.emplace_back(persistence_intervals);

            if (i == 1) {
//                persistence_intervals = {};
                for (int ip = 0; ip < already_paired_bigger_simplices.size(); ++ip) {
                    if (!already_paired_bigger_simplices[ip]) {
                        persistence_intervals.emplace_back(make_pair(
                                ((root->simplices[i - 1])[root->simplices[i - 1].size() - 1 - ip])->get_birth_time_protected(),
                                max_radii));
                    }
                }

            }
            persistence_diagram.emplace_back(persistence_intervals);


        }
        return persistence_diagram;
    }

    SubCloud(const Cloud &cloud, const vector<int> &points, float max_radius, int max_number_of_points_in_simplex) {
        max_radii = max_radius;
        matrix_of_distances = &cloud.distances;
        subcloud_of_points = points;
        root = new Simplex_tree(cloud.simplex_tree, cloud.simplex_tree->simplices.size());

        for (int i = 0; i < subcloud_of_points.size(); ++i) {
            //cout << i << " -> " << subcloud_of_points[i] << endl;
        }
        for (const auto &point: points) {
            auto found = root->find(vector<int>({point}), -1);
            root->simplices[0].emplace_back(found);
        }


        adjacency_matrix = vector<boost::dynamic_bitset<>>(points.size(), boost::dynamic_bitset<>(points.size()));
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
                if (cloud.distances[max(points[i], points[j])][min(points[i], points[j])] < max_radius) {
                    adjacency_matrix[i][j] = true;
                    if (points[i] < points[j]) {
//                        if (points[i] == 24 && points[j] == 27)
//                            exit(3);
                        root->insert(vector<int>({points[i], points[j]}),
                                     cloud.distances[max(points[i], points[j])][min(points[i], points[j])]);
                    }
                }
            }
        }
//        root->print();
    }

    friend ostream &operator<<(ostream &os, const vector<boost::dynamic_bitset<>> &matrix);
};


ostream &operator<<(ostream &os, const SubCloud &subcloud) {
    os << "Printing SubCloud" << endl;
    for (int i = 0; i < subcloud.adjacency_matrix.size(); ++i) {
        os << setw(space) << subcloud.subcloud_of_points[i] << ": ";
        for (int j = 0; j < subcloud.adjacency_matrix[i].size(); ++j) {
            os << subcloud.adjacency_matrix[i][j];
        }
        os << endl;
    }
    return os;
}


ostream &operator<<(ostream &os, const vector<vector<float>> &matrix) {
//    os << "Printing distances" << endl;
    for (const auto &row: matrix) {
        for (const auto &element: row) {
            os << element << " ";
        }
        os << endl;
    }
    return os;
}

template<typename T>
ostream &operator<<(ostream &os, const vector<T> &matrix) {
    for (const auto &row: matrix) {
        os << setw(space) << row;
    }
    //cout << endl;
    return os;
}





//template<typename T>
//ostream &operator<<(ostream &os, const vector<T> &matrix) {
////    os << "Printing distances" << endl;
//    for (const auto &row: matrix) {
//        os << setw(space) << row;
//    }
//    os << endl;
//    return os;
//}



vector<vector<pair<float, float>>> get_persistence_pairs(const Cloud &cloud, int count) {
    vector<int> dots;

    int max_number_of_points_in_simplex = 4;
    int number_of_dots = number_of_dots_in_file;
    set<int> not_include = {0, 2, 6, 7, 8, 9, 10, 23, 25};
    set<int> nnot_include = {0, 1, 2, 6, 8, 23, 24, 25};
    if (count == 0) {
        for (int i = 0; i < number_of_dots; ++i) {
            dots.emplace_back(i);
        }
    } else if (count == 1) {
        for (int i = 0; i < number_of_dots; ++i) {
            if (not_include.find(i) == not_include.end()) {
                dots.emplace_back(i);
            }
        }
    } else {
        for (int i = 0; i < number_of_dots; ++i) {
            if (nnot_include.find(i) == nnot_include.end()) {
                dots.emplace_back(i);
            }
        }
    }

    std::unique_ptr<SubCloud> subcloud(new SubCloud(cloud, dots, 8, max_number_of_points_in_simplex));
    mute.lock();
    std::cout << "started building tree" << std::endl;
    mute.unlock();
//
    subcloud->insert_all_simplices_including(max_number_of_points_in_simplex);
//    mute.lock();
//    std::cout << "inserted all simplices" << std::endl;
//    mute.unlock();
    vector<vector<pair<float,float>>> result;

    mute.lock();
    std::cout << "started calculating" << std::endl;
    mute.unlock();
    result = subcloud->calculate_betti_matrix(max_number_of_points_in_simplex - 1);
    mute.lock();
    for (int counter = result.size() - 1; counter >= 0; --counter) {
        cout << "dimension: " << result.size() - 1 - counter << endl;
        auto vector_of_pairs_of_dimension = result[counter];
        for (const auto &interval: vector_of_pairs_of_dimension) {
            cout << result.size() - 1 - counter << "         " << interval.first << ":" << interval.second << endl;
        }
    }
    for (int counter = result.size() - 1; counter >= 0; --counter) {
        cout << "dimension: " << result.size() - 1 - counter << " -< " << result[counter].size() << endl;
    }
    cout << "total errors " << number_of_errors << endl;
    mute.unlock();
    delete subcloud->root;
    return result;
}


int main() {
    //cout << "Hello, World!" << std::endl;


    int max_number_of_points_in_simplex = 4;
//    std::unique_ptr<Cloud> matrix(new Cloud("/home/leoneed/Desktop/dots1.off", max_number_of_points_in_simplex));
    std::unique_ptr<Cloud> matrix(new Cloud("/Users/leonardbee/CLionProjects/Persisitence_Diagram/dots_50", max_number_of_points_in_simplex));
//    std::unique_ptr<Cloud> matrix(new Cloud("/Users/leonardbee/Desktop/tore3D_60.off", max_number_of_points_in_simplex));


    std::ofstream fout("output.txt");
    fout << matrix->distances;
    vector<int> dots;

    int number_of_dots = number_of_dots_in_file;
    set<int> not_include = {0, 2, 6, 8, 23, 25};
    for (int i = 0; i < number_of_dots; ++i) {
        if (not_include.find(i) == not_include.end()) {
            dots.emplace_back(i);
        }
    }
    vector<std::thread> threads(5);

    for (int i = 0; i < threads.size(); ++i) {
        threads[i] = thread(get_persistence_pairs, *matrix, i);
    }
    for (int i = 0; i < threads.size(); ++i) {
        threads[i].join();
    }
    delete matrix->simplex_tree;




    return 0;
}
