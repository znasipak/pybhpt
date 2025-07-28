cmake_minimum_required(VERSION 3.15)
project(pybhpt LANGUAGES CXX)

# Find Python + NumPy
find_package(Python COMPONENTS Interpreter Development.Module NumPy REQUIRED)

# ---- Add options ----
set(PYBHPT_MARCH native CACHE STRING
                                 "Value of the -march compiler option if supported by compiler")

# Use Cython to compile the .pyx to .cpp
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

# Create the extension module
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
  cpp/src/unit_test.cpp
  WITH_SOABI
)

target_include_directories(cybhpt_full PRIVATE
  ${Python_NumPy_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include
)

set_property(TARGET cybhpt_full PROPERTY CXX_STANDARD 14)

# Compiler flags
include(CheckCXXCompilerFlag)
set(MARCH_OPT "-march=${PYBHPT_MARCH}")
check_cxx_compiler_flag("${MARCH_OPT}" COMPILER_SUPPORTS_MARCH)
if(COMPILER_SUPPORTS_MARCH)
  target_compile_options(cybhpt_full PRIVATE "${MARCH_OPT}")
endif()

# OpenMP (optional)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(cybhpt_full PRIVATE OpenMP::OpenMP_CXX)
endif()

# GSL + Boost
find_library(GSL_LIB gsl REQUIRED)
find_library(GSLCBLAS_LIB gslcblas REQUIRED)
find_library(BOOST_FILESYSTEM_LIB boost_filesystem REQUIRED)
target_link_libraries(cybhpt_full PRIVATE ${GSL_LIB} ${GSLCBLAS_LIB} ${BOOST_FILESYSTEM_LIB})

# Install (optional)
install(TARGETS cybhpt_full DESTINATION "./")
