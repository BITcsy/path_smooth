cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(Qp-smooth)

# dependences & Library target
find_package(OsqpEigen)
find_package(Eigen3)
find_package(PythonLibs 3.6)
 
# Examples
set(SRC_LIST 
    qp_smooth.cpp 
    scenarios.cpp)

include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR} "/home/csy/.local/lib/python3.6/site-packages/numpy/core/include")
# include_directories("/mnt/d/linux_workspace/csy1_workspace/src" "/home/csy/.local/lib/python3.6/site-packages/numpy/core/include")
add_executable(qp_smooth ${SRC_LIST})
target_include_directories(qp_smooth PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(qp_smooth ${PYTHON_LIBRARIES} OsqpEigen::OsqpEigen)

# target_link_libraries(minimal PRIVATE matplotlib_cpp)
# set_target_properties(minimal PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")


