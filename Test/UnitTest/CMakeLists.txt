cmake_minimum_required(VERSION 3.10)

project(UnitTest)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(MATRIX_INCLUDE_DIR ../../include)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")

add_executable(${PROJECT_NAME} test.cpp)

target_include_directories(${PROJECT_NAME} PRIVATE ${MATRIX_INCLUDE_DIR})

target_link_libraries(${PROJECT_NAME} gtest)