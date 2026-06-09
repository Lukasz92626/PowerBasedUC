cmake_minimum_required(VERSION 3.8)
project(PowerBasedUC)
set(CMAKE_CXX_STANDARD 20)

set(CPLEX_ROOT_PATH "C:/Program Files/IBM/ILOG/CPLEX_Studio222")

add_executable(PowerBasedUC "PowerBasedUC.cpp")

target_include_directories(PowerBasedUC PUBLIC
    "${CPLEX_ROOT_PATH}/concert/include"
    "${CPLEX_ROOT_PATH}/cplex/include"
)
target_link_directories(PowerBasedUC PUBLIC
    "${CPLEX_ROOT_PATH}/concert/lib/x64_windows_msvc14/stat_mdd"
    "${CPLEX_ROOT_PATH}/cplex/lib/x64_windows_msvc14/stat_mdd"
)
target_link_libraries(PowerBasedUC PUBLIC concert cplex2220 ilocplex)

if(WIN32)
    target_compile_definitions(PowerBasedUC PRIVATE _CRT_SECURE_NO_WARNINGS)
endif()