add_executable(eskf_example eskf_example.cpp ESKF.cpp math_utils.cpp)
target_link_libraries(eskf_example DataLoader utils Eigen3::Eigen ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})