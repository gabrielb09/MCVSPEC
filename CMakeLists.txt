cmake_minimum_required(VERSION 3.15)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(mcvspec VERSION 2.0.0)


cmake_minimum_required(VERSION 3.16)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(MCVSPEC LANGUAGES CXX VERSION 2.0.0)

add_library(mcvspec STATIC source/Cataclysmic_Variable.cpp)

target_compile_features(mcvspec PUBLIC cxx_std_20)

target_include_directories(mcvspec PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include/mcvspec
  ${CMAKE_CURRENT_SOURCE_DIR}/include/mcvspec/integrators
  ${CMAKE_CURRENT_SOURCE_DIR}/include/mcvspec/root_finders
  ${CMAKE_CURRENT_SOURCE_DIR}/resources
)

# Tests
enable_testing()

add_executable(test_rk45 tests/rk45_test.cpp)
add_executable(test_mcv tests/mcv_test.cpp)

target_link_libraries(test_rk45 PRIVATE mcvspec)
target_link_libraries(test_mcv PRIVATE mcvspec)

add_test(NAME rk45 COMMAND test_rk45)
add_test(NAME mcv COMMAND test_mcv)
