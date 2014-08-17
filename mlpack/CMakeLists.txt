project(mlpack)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    "${CMAKE_SOURCE_DIR}/cmake/Modules")

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIR})
set(LIBS ${LIBS} ${ARMADILLO_LIBRARIES})

find_package(Mlpack REQUIRED)
include_directories(${Mlpack_INCLUDE_DIRS})
set(LIBS ${LIBS} ${Mlpack_LIBRARIES})

find_package(LibXml2 REQUIRED)
include_directories(${LIBXML2_INCLUDE_DIR})
set(LIBS ${LIBS} ${LIBXML2_LIBRARIES})

add_definitions(-std=c++11)
add_definitions(-DARMA_USE_CXX11)

aux_source_directory(. SRC_LIST)

add_executable(${PROJECT_NAME} ${SRC_LIST})
target_link_libraries(${PROJECT_NAME} ${LIBS})

