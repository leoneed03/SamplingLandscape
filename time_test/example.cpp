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
    /// 3. double - max radii of VR complex
    /// 4. int - num of threads >= 1
    /// 5. int - num of samples >= 1
    /// 6. density coefficient (0, 1]

    double time1 = landscape_ripser(path, "results", 2, radii, 1, 1, 1);
    double time2 = landscape_algorithm(path, "landscapes", 2, radii, 1, 1, 0.9);
    double time3 = landscape_gudhi(path, "landscapes", 2, radii, 1, 1, 1);

    return 0;
}