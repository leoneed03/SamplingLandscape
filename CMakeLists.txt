cmake_minimum_required(VERSION 3.10)
project(restribution)

set(CMAKE_CXX_STANDARD 17)

FIND_PACKAGE(Boost)
MESSAGE("Boost information:")
MESSAGE("  Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
MESSAGE("  Boost_LIBRARIES: ${Boost_LIBRARIES}")
MESSAGE("  Boost_LIBRARY_DIRS: ${Boost_LIBRARY_DIRS}")

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})


list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})


find_package(TBB)
include_directories(${TBB_INCLUDE_DIRS})
add_definitions(${TBB_DEFINITIONS})


LINK_DIRECTORIES(${tbb_LIBRARY_DIRS})

add_executable(restribution main.cpp)
target_link_libraries(restribution ${TBB_LIBRARIES})