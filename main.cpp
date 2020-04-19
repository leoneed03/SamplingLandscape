#include <iostream>
#include <random>
#include <iterator>
#include "ripser.cpp"

using namespace std;

set<int> get_random_sample(vector<int>& vector_of_points, int size_of_one_sample) {
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(vector_of_points.begin(), vector_of_points.end(), g);
    for (const auto& e: vector_of_points) {
        cout << e << " ";
    }
    cout << endl;
    sort(vector_of_points.begin(), vector_of_points.begin() + size_of_one_sample);
//    ripser_func();
    std::set<int> subcloud(vector_of_points.begin(), vector_of_points.begin() + size_of_one_sample);

    return std::set<int>(vector_of_points.begin(), vector_of_points.begin() + size_of_one_sample);
}


int main() {
    vector<int> v;
    int n = 10;

    for (int i = 0; i < n; ++i) {
//        v.emplace_back(i);
    }
    for (int i = 0; i < n; ++i) {
//        get_random_sample(v);
    }


    double max_edge_length = 8;
    double min_edge_length = 1e-10;
    int max_cohomology_rank = 2;
    std::string filename = "/Users/leonardbee/CLionProjects/restribution/cloud_50";
    std::ifstream file_stream(filename);
    int number_of_points = 0;
    std::string line;
    while (std::getline(file_stream, line)) {
        ++number_of_points;
        std::vector<value_t> point;
        std::istringstream s(line);
    }
    std::cout << "total number of points " << number_of_points << std::endl;

    vector<int> cloud(number_of_points);
    for (int i = 0; i < number_of_points; ++i) {
        cloud[i] = i;
    }
    std::vector<std::string> argv_strings = {"./ripser",  "--dim",  "2",  "--threshold",  "8",
                                             "--modulus",  "2", "--format", "point-cloud", "/Users/leonardbee/CLionProjects/restribution/cloud_50"};
    auto res = main_ripser(10, argv_strings, get_random_sample(cloud, (int) (cloud.size() * 1)));
    std::cout << res.size();
    for (const auto& dim: res) {
        for (const auto& persistence_pair: dim) {
            std::cout << persistence_pair.first << " vs " << persistence_pair.second << std::endl;
        }
        std::cout << "____________________________________\n" << std::endl;
    }

    return 0;
}