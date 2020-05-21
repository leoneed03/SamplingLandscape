# Subsampling Methods For Persistence Landscape

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
or by commands 
```
sudo cp -r thread_pool/thread_pool /usr/local/include/
sudo cp -r thread_pool/thread_pool.hpp /usr/local/include/
```

## Usage
see time_test/example.cpp

## Dependencies

1. GUDHI: https://github.com/GUDHI/gudhi-devel
2. PHAT: https://github.com/blazs/phat
3. TBB: https://github.com/oneapi-src/oneTBB
4. Boost: https://www.boost.org/users/download/
5. Thread pool: https://github.com/inkooboo/thread-pool-cpp
