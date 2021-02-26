cmake_minimum_required(VERSION 3.1)
set (CMAKE_CXX_STANDARD 11)
project(Region3D)
find_package(CGAL)
include( ${CGAL_USE_FILE} )
add_executable(Region3D UF.h main.cpp)
add_executable(Refine refine.cpp)


