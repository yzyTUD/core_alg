cmake_minimum_required(VERSION 2.6)
project("zrlib")

include_directories(../zrlib_cpp_headers/)
add_executable(test test.cpp)
configure_file(1.txt 1.txt COPYONLY)