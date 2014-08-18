project(custom)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)

add_definitions(-std=c++11)

#find_package(
#    Boost REQUIRED
#    Components multiarray
#)

add_executable(${PROJECT_NAME} ${SRC_LIST})

