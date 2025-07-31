cmake_minimum_required(VERSION 3.15)
project(pybhpt LANGUAGES CXX)

# Python & NumPy
find_package(Python COMPONENTS Interpreter Development.Module NumPy REQUIRED)

# Conda detection
if(DEFINED ENV{CONDA_PREFIX})
  message(STATUS "Using Conda environment at $ENV{CONDA_PREFIX}")
  list(APPEND CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}")
  list(APPEND CMAKE_LIBRARY_PATH "$ENV{CONDA_PREFIX}/lib")
  list(APPEND CMAKE_INCLUDE_PATH "$ENV{CONDA_PREFIX}/include")
endif()

# ------------------------- GSL Detection -------------------------
find_package(GSL REQUIRED)

# ------------------------- Boost Headers Only -------------------------
set(BOOST_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/extern/boost")

# ------------------------- Cython -------------------------
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/cython")

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/cython/redshift_wrap.cpp"
  COMMAND ${Python_EXECUTABLE} -m cython
          "${CMAKE_CURRENT_SOURCE_DIR}/cython/redshift_wrap.pyx"
          --output-file "${CMAKE_CURRENT_BINARY_DIR}/cython/redshift_wrap.cpp"
          -3 -+ --module-name "cybhpt_full"
          -I "${CMAKE_CURRENT_SOURCE_DIR}/cpp/include"
  DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/cython/redshift_wrap.pyx"
  COMMENT "Generating redshift_wrap.cpp from Cython"
  VERBATIM
)

# ------------------------- Extension Module -------------------------
python_add_library(cybhpt_full MODULE
  "${CMAKE_CURRENT_BINARY_DIR}/cython/redshift_wrap.cpp"
  cpp/src/utils.cpp
  cpp/src/specialfunc.cpp
  cpp/src/kerr.cpp
  cpp/src/geo.cpp
  cpp/src/cf.cpp
  cpp/src/swsh.cpp
  cpp/src/monodromy.cpp
  cpp/src/hypergeo_f.cpp
  cpp/src/bessel.cpp
  cpp/src/hypergeo_u.cpp
  cpp/src/nusolver.cpp
  cpp/src/mst.cpp
  cpp/src/gsn_asymp.cpp
  cpp/src/radialsolver.cpp
  cpp/src/sourceintegration.cpp
  cpp/src/teukolsky.cpp
  cpp/src/hertz.cpp
  cpp/src/metriccoeffs.cpp
  cpp/src/fluxes.cpp
  cpp/src/metric.cpp
  cpp/src/redshift.cpp
  WITH_SOABI
)

# ------------------------- Include Paths -------------------------
target_include_directories(cybhpt_full PRIVATE
  ${Python_NumPy_INCLUDE_DIRS}
  ${BOOST_INCLUDE_DIR}
  ${GSL_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include
)

# C++ flags
set_property(TARGET cybhpt_full PROPERTY CXX_STANDARD 14)
set(PYBHPT_MARCH native CACHE STRING "Value of the -march compiler option")
include(CheckCXXCompilerFlag)
set(MARCH_OPT "-march=${PYBHPT_MARCH}")
check_cxx_compiler_flag("${MARCH_OPT}" COMPILER_SUPPORTS_MARCH)
if(COMPILER_SUPPORTS_MARCH)
  target_compile_options(cybhpt_full PRIVATE "${MARCH_OPT}")
endif()

# OpenMP
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(cybhpt_full PRIVATE OpenMP::OpenMP_CXX)
else()
  message(STATUS "OpenMP not found, proceeding without it")
  target_compile_definitions(cybhpt_full PRIVATE PYBHPT_NO_OPENMP)
endif()

# GSL
target_link_libraries(cybhpt_full PRIVATE GSL::gsl GSL::gslcblas)
set_target_properties(cybhpt_full PROPERTIES
  INSTALL_RPATH "$ORIGIN/../../opt/gsl-install/lib"
  BUILD_WITH_INSTALL_RPATH TRUE
)


# Optional install
install(TARGETS cybhpt_full DESTINATION "./")
