cmake_minimum_required(VERSION 3.12)
project(cppmodel)

find_package(OpenMP REQUIRED)
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
find_package(BLAS)
message( STATUS BLAS found: ${BLAS_LIBRARIES} )
find_package(LAPACK)
message( STATUS LAPACK found: ${LAPACK_LIBRARIES} )

set(CMAKE_CXX_STANDARD 11)

add_definitions(-DSTATS_USE_ARMA)
add_definitions(-DARMA_DONT_USE_WRAPPER)
#add_definitions(-DARMA_EXTRA_DEBUG)


add_library(cppmodel SHARED
            blocked_gibbs_imp.hpp blocked_gibbs.hpp
            cluster.hpp options.hpp
            semi_hdp_hierarchy.hpp semi_hdp_hierarchy.cpp
            random_engine.hpp random_engine.cpp
            utils.hpp utils.cpp
            statslib/stats.hpp
            metropolis_hastings.hpp
            marginal_semi_hdp.hpp marginal_semi_hdp_imp.hpp
            memory_collector_by_param.hpp memory_collector_by_param.cpp
            memory_collector_by_iter.hpp memory_collector_by_iter.cpp)

set_target_properties(cppmodel PROPERTIES VERSION 0.1)
target_include_directories(cppmodel PRIVATE .)
target_link_libraries(cppmodel
      ${ARMADILLO_LIBRARIES}
      ${BLAS_LIBRARIES}
      ${LAPACK_LIBRARIES})

install(TARGETS cppmodel DESTINATION ./)

FILE(GLOB CppSources RELATIVE ${CMAKE_SOURCE_DIR} "*.cpp")
list(REMOVE_ITEM CppSources rcpp_functions.cpp)
list(REMOVE_ITEM CppSources RcppExports.cpp)


add_executable(test_blocked_gibbs
        spikes/test_blocked_gibbs.cpp
        utils.cpp
        ${CppSources})

target_link_libraries(test_blocked_gibbs
        PUBLIC OpenMP::OpenMP_CXX
        ${ARMADILLO_LIBRARIES}
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES})

add_executable(test_mh_update
        spikes/test_mh_update.cpp
        utils.cpp
        semi_hdp_hierarchy.cpp
        semi_hdp_hierarchy.hpp
        ${CppSources})

add_executable(test_stirling
        spikes/test_stirling.cpp
        utils.cpp)

target_link_libraries(test_mh_update
        PUBLIC OpenMP::OpenMP_CXX
        ${ARMADILLO_LIBRARIES}
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES})

add_executable(test_conditional
        spikes/test_conditional_semihdp.cpp
        utils.cpp
        ${CppSources})


target_link_libraries(test_conditional
        PUBLIC OpenMP::OpenMP_CXX
        ${ARMADILLO_LIBRARIES}
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES})

add_executable(test_nig
        spikes/test_nig_update.cpp
        utils.cpp
        ${CppSources})

target_link_libraries(test_nig
        PUBLIC OpenMP::OpenMP_CXX
        ${ARMADILLO_LIBRARIES}
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES})
