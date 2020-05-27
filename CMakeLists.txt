set(CMAKE_C_COMPILER gcc)
cmake_minimum_required(VERSION 3.1)

project(sampling)

include (CTest)
set(CMAKE_CXX_STANDARD 14)
set(GCC_COVERAGE_COMPILE_FLAGS " -Os -DNDEBUG")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
#set(CMAKE_EXE_LINKER_FLAGS "-lpthread -lboost_system ")

FIND_PACKAGE(Boost COMPONENTS thread system filesystem REQUIRED)
MESSAGE("Boost information:")
MESSAGE("   Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
MESSAGE("   Boost_LIBRARIES: ${Boost_LIBRARIES}")
MESSAGE("   Boost_LIBRARY_DIRS: ${Boost_LIBRARY_DIRS}")

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
include_directories(/usr/local/include)
find_package(Threads)

LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

find_package(TBB)
include_directories(${TBB_INCLUDE_DIRS})
add_definitions(${TBB_DEFINITIONS})


LINK_DIRECTORIES(${tbb_LIBRARY_DIRS})


find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_executable(sampling_algorithm time_test/test_algorithm.cpp)
target_link_libraries(sampling_algorithm ${GTEST_LIBRARIES} ${TBB_LIBRARIES} ${Boost_LIBRARIES} pthread)
add_test(NAME sampling_algorithm
        COMMAND sampling_algorithm)

add_executable(sampling_magn test/test_diagram_computation.cpp)
target_link_libraries(sampling_magn ${GTEST_LIBRARIES} ${TBB_LIBRARIES} ${Boost_LIBRARIES} pthread)
add_test(NAME sampling_magn
        COMMAND sampling_magn)

add_executable(sampling_1 test/test_diagram_computation_v.cpp)
target_link_libraries(sampling_1 ${GTEST_LIBRARIES} ${TBB_LIBRARIES} ${Boost_LIBRARIES} pthread)
add_test(NAME sampling_1
        COMMAND sampling_1)

add_executable(sampling_gudhi time_test/test_gudhi.cpp)
target_link_libraries(sampling_gudhi ${GTEST_LIBRARIES} ${TBB_LIBRARIES} ${Boost_LIBRARIES} pthread)
add_test(NAME sampling_gudhi
        COMMAND sampling_gudhi)



add_executable(sampling_ripser time_test/test_ripser.cpp)
target_link_libraries(sampling_ripser ${GTEST_LIBRARIES} ${TBB_LIBRARIES} ${Boost_LIBRARIES} pthread)
add_test(NAME sampling_ripser
        COMMAND sampling_ripser)
