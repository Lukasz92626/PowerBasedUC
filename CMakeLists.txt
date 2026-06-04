cmake_minimum_required(VERSION 3.8)

project(PowerBasedUC)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Ścieżka do CPLEX (dostosuj do swojej wersji)
set(CPLEX_ROOT_PATH "C:/Program Files/IBM/ILOG/CPLEX_Studio_Community2212" CACHE PATH "Path to CPLEX installation")

# Automatycznie dodaj wszystkie pliki .cpp i .h
file(GLOB SOURCES "*.cpp")
file(GLOB HEADERS "*.h")

# Dodaj plik wykonywalny
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS} "PowerBasedUC.cpp")

# Dołącz nagłówki CPLEX
target_include_directories(${PROJECT_NAME} PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CPLEX_ROOT_PATH}/concert/include"
    "${CPLEX_ROOT_PATH}/cplex/include"
)

# Dołącz biblioteki CPLEX
target_link_directories(${PROJECT_NAME} PRIVATE
    "${CPLEX_ROOT_PATH}/concert/lib/x64_windows_msvc14/stat_mdd"
    "${CPLEX_ROOT_PATH}/cplex/lib/x64_windows_msvc14/stat_mdd"
)

# Linkuj biblioteki
target_link_libraries(${PROJECT_NAME}
    concert
    cplex2212
    ilocplex
)

# Dla Windows
if(WIN32)
    target_compile_definitions(${PROJECT_NAME} PRIVATE _CRT_SECURE_NO_WARNINGS)
endif()