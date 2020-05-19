# SubsamplingMethodsForPersistenceLandscape

## Other sources (Linux)

```
sudo apt-get -y install libtbb-dev
sudo apt-get -y install libgudhi-dev
sudo apt-get -y install libphat-dev
```
Thread pool can be installed from https://github.com/inkooboo/thread-pool-cpp
or by commands 
```
- sudo cp -r thread_pool/thread_pool /usr/local/include/
- sudo cp -r thread_pool/thread_pool.hpp /usr/local/include/
```

## Usage

```c++
cd time_test
g++ test_{version}.cpp -std=c++17 -ltbb [-Xpreprocessor] -fopenmp -Os -lomp -o algo
./algo
```
## Dependencies

1. GUDHI: https://github.com/GUDHI/gudhi-devel
2. PHAT: https://github.com/blazs/phat
3. TBB: https://github.com/oneapi-src/oneTBB
4. Boost: https://www.boost.org/users/download/
