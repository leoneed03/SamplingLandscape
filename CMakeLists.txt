set(CMAKE_C_COMPILER gcc)
cmake_minimum_required (VERSION 3.1)

project (sampling)

set(CMAKE_CXX_STANDARD 14)
set(GCC_COVERAGE_COMPILE_FLAGS " -Os -DNDEBUG")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )


INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
include_directories(/usr/local/include)

LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

find_package(TBB)
include_directories(${TBB_INCLUDE_DIRS})
add_definitions(${TBB_DEFINITIONS})


LINK_DIRECTORIES(${tbb_LIBRARY_DIRS})


find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_executable(sampling_1 test/test_diagram_computation_v.cpp)
#add_executable(sampling test/test_diagram_computation_v_1.cpp)
target_link_libraries(sampling_1 ${GTEST_LIBRARIES} ${TBB_LIBRARIES} pthread)
#target_link_libraries(sampling_1 /usr/local/lib/libgtest.a ${TBB_LIBRARIES})

enable_testing()
add_test(NAME    sampling_1
        COMMAND sampling_1)

