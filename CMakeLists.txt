cmake_minimum_required(VERSION 3.15)
project(pybhpt LANGUAGES CXX)

# Find Python + NumPy
find_package(Python COMPONENTS Interpreter Development.Module NumPy REQUIRED)

if(DEFINED ENV{CONDA_PREFIX})
  message(STATUS "Using Conda environment at $ENV{CONDA_PREFIX}")
  list(APPEND CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}")
  list(APPEND CMAKE_LIBRARY_PATH "$ENV{CONDA_PREFIX}/lib")
  list(APPEND CMAKE_INCLUDE_PATH "$ENV{CONDA_PREFIX}/include")
endif()

# ------------------------- Boost + GSL Setup -------------------------

if(DEFINED CMAKE_TOOLCHAIN_FILE)
  message(STATUS "Using vcpkg toolchain: ${CMAKE_TOOLCHAIN_FILE}")
  
  # Use vcpkg to find Boost + GSL
  find_package(GSL REQUIRED)
  if(GSL_FOUND)
    message(STATUS "GSL version: ${GSL_VERSION}")
  else()
    message(FATAL_ERROR "GSL not found via vcpkg")
  endif()

  find_package(Boost REQUIRED COMPONENTS headers)
  if(Boost_FOUND)
    message(STATUS "Boost version: ${Boost_VERSION}")
  else()
    message(FATAL_ERROR "Boost not found via vcpkg")
  endif()

else()
  message(WARNING "vcpkg toolchain not detected. Falling back to manual linking.")

  # Manual fallback if not using vcpkg
  find_library(GSL_LIB gsl REQUIRED)
  find_library(GSLCBLAS_LIB gslcblas REQUIRED)
  find_library(BOOST_FILESYSTEM_LIB boost_filesystem REQUIRED)

  if(NOT GSL_LIB OR NOT GSLCBLAS_LIB)
    message(FATAL_ERROR "Could not find GSL or GSLCBLAS libraries")
  else()
    message(STATUS "Found GSL version: ${GSL_LIB}")
    message(STATUS "Found GSLCBLAS version: ${GSLCBLAS_LIB}")
  endif()

  if(NOT BOOST_FILESYSTEM_LIB)
    message(FATAL_ERROR "Could not find Boost Filesystem library")
  else()
    message(STATUS "Found Boost Filesystem version: ${BOOST_FILESYSTEM_LIB}")
  endif()
endif()

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
else()
  message(STATUS "OpenMP not found, proceeding without it")
  # Optionally, you can define a macro to indicate OpenMP is not used
  target_compile_definitions(cybhpt_full PRIVATE PYBHPT_NO_OPENMP)
endif()

if(DEFINED CMAKE_TOOLCHAIN_FILE)
  target_link_libraries(cybhpt_full PRIVATE GSL::gsl GSL::gslcblas)
  target_link_libraries(cybhpt_full PRIVATE Boost::headers)
  get_target_property(boost_includes Boost::headers INTERFACE_INCLUDE_DIRECTORIES)
message(STATUS "Boost::headers includes: ${boost_includes}")
else()
  target_link_libraries(cybhpt_full PRIVATE ${GSL_LIB} ${GSLCBLAS_LIB} ${BOOST_FILESYSTEM_LIB})
endif()


# Install (optional)
install(TARGETS cybhpt_full DESTINATION "./")
