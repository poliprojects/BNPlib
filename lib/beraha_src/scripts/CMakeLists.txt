cmake_minimum_required(VERSION 3.12)

project(experiments)

find_package(OpenMP)
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
find_package(BLAS)
message( STATUS BLAS found: ${BLAS_LIBRARIES} )
find_package(LAPACK)
message( STATUS LAPACK found: ${LAPACK_LIBRARIES} )

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffast-math -funsafe-math-optimizations -ftracer")

add_definitions(-DSTATS_USE_ARMA)
add_definitions(-DARMA_DONT_USE_WRAPPER)
#add_definitions(-DARMA_EXTRA_DEBUG)


FILE(GLOB CppSources RELATIVE ${CMAKE_SOURCE_DIR} "../src/*.cpp")
message(CppSources: ${CppSources})
list(REMOVE_ITEM CppSources ../src/rcpp_functions.cpp)
list(REMOVE_ITEM CppSources ../src/RcppExports.cpp)


add_executable(two_population_experiments
        two_population_experiments.cpp
        ../src/utils.cpp
        ${CppSources})

add_executable(two_population_meta
        two_pouplation_meta_experiment.cpp
        ../src/utils.cpp
        ${CppSources})

add_executable(performance_comparison
        compare_performance.cpp
        ../src/utils.cpp
        ${CppSources})

add_executable(four_population
  four_population_experiment.cpp
  ../src/utils.cpp
  ${CppSources}
)

add_executable(many_population
        many_population_experiment.cpp
        ../src/utils.cpp
        ${CppSources})

add_executable(grades_chile
  grades_chile.cpp
  ../src/utils.cpp
  ${CppSources}
)


if(OpenMP_CXX_FOUND)
    target_link_libraries(two_population_experiments
      PUBLIC OpenMP::OpenMP_CXX
      ${ARMADILLO_LIBRARIES}
      ${BLAS_LIBRARIES}
      ${LAPACK_LIBRARIES})

    target_link_libraries(two_population_meta
      PUBLIC OpenMP::OpenMP_CXX
      ${ARMADILLO_LIBRARIES}
      ${BLAS_LIBRARIES}
      ${LAPACK_LIBRARIES})

    target_link_libraries(performance_comparison
      PUBLIC OpenMP::OpenMP_CXX
      ${ARMADILLO_LIBRARIES}
      ${BLAS_LIBRARIES}
      ${LAPACK_LIBRARIES})

    target_link_libraries(four_population
      PUBLIC OpenMP::OpenMP_CXX
      ${ARMADILLO_LIBRARIES}
      ${BLAS_LIBRARIES}
      ${LAPACK_LIBRARIES})

    target_link_libraries(many_population
      PUBLIC OpenMP::OpenMP_CXX
      ${ARMADILLO_LIBRARIES}
      ${BLAS_LIBRARIES}
      ${LAPACK_LIBRARIES})

    target_link_libraries(grades_chile
      PUBLIC OpenMP::OpenMP_CXX
      ${ARMADILLO_LIBRARIES}
      ${BLAS_LIBRARIES}
      ${LAPACK_LIBRARIES})
else()
    target_link_libraries(two_population_experiments
        ${ARMADILLO_LIBRARIES}
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES})

    target_link_libraries(two_population_meta
        ${ARMADILLO_LIBRARIES}
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES})

    target_link_libraries(performance_comparison
      ${ARMADILLO_LIBRARIES}
      ${BLAS_LIBRARIES}
      ${LAPACK_LIBRARIES})

    target_link_libraries(four_population
      ${ARMADILLO_LIBRARIES}
      ${BLAS_LIBRARIES}
      ${LAPACK_LIBRARIES})

    target_link_libraries(many_population
      ${ARMADILLO_LIBRARIES}
      ${BLAS_LIBRARIES}
      ${LAPACK_LIBRARIES})

    target_link_libraries(grades_chile
      ${ARMADILLO_LIBRARIES}
      ${BLAS_LIBRARIES}
      ${LAPACK_LIBRARIES})
endif()
