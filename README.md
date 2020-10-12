# Subsampling Methods For Vietoris-Rips Persistence Landscape

## Libraries installation (Linux)

TBB, GUDHI, PHAT can be installed with package manager
```
sudo apt-get -y install libtbb-dev
sudo apt-get -y install libgudhi-dev
sudo apt-get -y install libphat-dev
```
or directly from their repositories (see "Dependencies") 
or manually from "dep"

Thread pool can be installed from https://github.com/inkooboo/thread-pool-cpp

## Usage

every implemented landscape-building algorithm has same interface
let's use landscape_ripser from ripser/landscape_r.h:

Function arguments should be:
  std::string from -- path to file point cloud descriptor where each line represents vector of coordinates for one point (see dataset folder for examples)
  std::string to -- path to file where to store resulting landscape
  int max_rank -- max rank of homology computated
  double max_edge_length -- threshold for edge in connectivity graph 
  int number_of_thread_workers = 1 -- thread pool size
  int number_of_samples = 1 -- number of samples (more for more accuracy)
  double subsample_density_coefficient = 1.0 -- fraction of sample size and cloud size (less for less computation time)

## Example

for each launch of landscape_{method} you will get set of files with landscapes which can be easily represented:

![Image of Yaktocat](https://drive.google.com/file/d/1EIh2edGS-uiTFr45ZFREAlvYyg1lUex_/view?usp=sharing)

left -- exact landscape, right -- approximated one with subsample_density_coefficient = 0.4 and number_of_samples = 10
## Dependencies

1. GUDHI: https://github.com/GUDHI/gudhi-devel
2. PHAT: https://github.com/blazs/phat
3. TBB: https://github.com/oneapi-src/oneTBB
4. Boost: https://www.boost.org/users/download/
5. Thread pool: https://github.com/inkooboo/thread-pool-cpp
