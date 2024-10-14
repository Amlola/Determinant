cmake_minimum_required(VERSION 3.10)

project(EndToEndTest)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(MATRIX_INCLUDE_DIR ../../include)

add_executable(${PROJECT_NAME} main.cpp)

target_include_directories(${PROJECT_NAME} PRIVATE ${MATRIX_INCLUDE_DIR})