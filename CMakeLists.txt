cmake_minimum_required(VERSION 3.18)                                                                               
project(powersasa_test)

add_executable(sasatest
sasatest.cpp
power_diagram.h
power_sasa.h
)

add_definitions(-D_GLIBCXX_NO_ASSERTIONS)
set_property(TARGET sasatest PROPERTY CXX_STANDARD 20)

include_directories(".")
