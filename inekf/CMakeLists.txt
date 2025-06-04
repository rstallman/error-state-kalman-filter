cmake_minimum_required(VERSION 3.10)
project(modrob CXX)

include(CMakeHeader.cmake)

find_package(Eigen3 REQUIRED)

find_package(Boost REQUIRED COMPONENTS unit_test_framework)

# include_directories(${EIGEN3_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

# ##############################################################################
# IEKF
add_library(iekf SHARED src/IEKF.cpp)
target_link_libraries(iekf utils Eigen3::Eigen ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
target_add_test_case(test_iekf test_iekf.cpp)

target_link_libraries(test_iekf iekf Eigen3::Eigen ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

# ################################################################################
# ESKF
# add_library(eskf SHARED src/ESKF.cpp)
# target_link_libraries(eskf utils Eigen3::Eigen ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

# ##############################################################################

# ##############################################################################
# DataLoader
add_library(DataLoader SHARED src/DataLoader.cpp)
add_executable(iekf_example src/iekf_example.cpp)
target_link_libraries(iekf_example DataLoader iekf Eigen3::Eigen ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

# add_executable(eskf_example src/eskf_example.cpp)
# target_link_libraries(eskf_example DataLoader eskf Eigen3::Eigen ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_executable(eskf_deepseek src/eskf_deepseek.cpp)
target_link_libraries(eskf_deepseek DataLoader utils Eigen3::Eigen ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

# ##############################################################################

# ##############################################################################
# utils
add_library(utils SHARED src/utils.cpp)
target_add_test_case(test_utils test_utils.cpp)
target_link_libraries(test_utils utils Eigen3::Eigen ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

# ##############################################################################
