#include "../include/algorithm/landscape_a.h"
#include "../include/ripser/landscape_r.h"
#include "../include/gudhi/landscape_g.h"
#include "../util/compare.cpp"
#include <iostream>


using namespace smpl;

int main() {
    /// path to file with dots
    /// file format: one line - one point
    /// coordinates in line separated with spaces



    std::string path = "../dataset/magnetometer/s50.txt";
    double radii = 1e5;

    /// args:
    /// 1. std::string - path to file with points
    /// 2. std::string - path where to save landscape
    /// 3. max_diagram_rank
    /// 4. double - max radii of VR complex
    /// 5. int - num of threads >= 1
    /// 6. int - num of samples >= 1
    /// 7. density coefficient (0, 1]


    tbb::concurrent_vector<tbb::concurrent_vector<std::vector<std::pair<double, double>>>> apd;

    double time1 = landscape_ripser(path, "results", 2, radii, 3, 9, 1);
    double time2 = landscape_algorithm(path, "landscapes", 2, radii, 1, 2, 0.9);
    double time3 = landscape_gudhi(path, "landscapes", 2, radii, 4, 10, 0.5);
//    double time1 = landscape_ripser(path, "results", 2, radii, 4, 5, 0.9);
//    double time3 = landscape_algorithm_with_diagrams(path, "landscapes", apd, 2, radii, 4, 15, 0.9);
//    double time2 = landscape_algorithm(path, "landscapes", 2, radii, 4, 15, 0.9);

//    std::cout << "S_amples " << apd.size() << std::endl;
//    for (const auto& e: apd) {
//        std::cout << "            dims " << e.size() << std::endl;
//        for (const auto& a: e) {
//            std::cout << "                         intervals " << a.size() << std::endl;
//        }
//    }


    return 0;
}