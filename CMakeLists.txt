cmake_minimum_required(VERSION 3.10)
project(innovate_fpga)
find_package( OpenCV REQUIRED )
set(CMAKE_CXX_STANDARD 11)

add_executable(innovate_fpga main.cpp)

target_link_libraries( innovate_fpga ${OpenCV_LIBS} )
