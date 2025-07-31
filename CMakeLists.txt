cmake_minimum_required(VERSION 3.15)
project(pybhpt LANGUAGES CXX)

# Find Python + NumPy
find_package(Python COMPONENTS Interpreter Development.Module NumPy REQUIRED)

# Auto-detect Conda and set search paths
if(DEFINED ENV{CONDA_PREFIX})
  message(STATUS "Using Conda environment at $ENV{CONDA_PREFIX}")
  list(APPEND CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}")
  list(APPEND CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}/lib/cmake")
  list(APPEND CMAKE_LIBRARY_PATH "$ENV{CONDA_PREFIX}/lib")
  list(APPEND CMAKE_INCLUDE_PATH "$ENV{CONDA_PREFIX}/include")
endif()

# ------------------------- Boost + GSL Setup -------------------------

if(DEFINED CMAKE_TOOLCHAIN_FILE)
  message(STATUS "Using vcpkg toolchain: ${CMAKE_TOOLCHAIN_FILE}")
  find_package(GSL REQUIRED)
  if(GSL_FOUND)
    message(STATUS "GSL version: ${GSL_VERSION}")
  else()
    message(FATAL_ERROR "GSL not found via vcpkg")
  endif()
  # Assume Boost was installed and found by vcpkg
  find_package(Boost REQUIRED COMPONENTS math numeric_odeint)
  if(Boost_FOUND)
    message(STATUS "Boost version: ${Boost_VERSION}")
  else()
    message(FATAL_ERROR "Boost not found via vcpkg")
  endif()
else()
  message(WARNING "vcpkg toolchain not detected. Falling back to FetchContent Boost and manual GSL")

  include(FetchContent)

  # --------- Find GSL manually ---------
  find_path(GSL_INCLUDE_DIR gsl/gsl_math.h PATHS ${CMAKE_INCLUDE_PATH} REQUIRED)
  find_library(GSL_LIB NAMES gsl PATHS ${CMAKE_LIBRARY_PATH} REQUIRED)
  find_library(GSLCBLAS_LIB NAMES gslcblas PATHS ${CMAKE_LIBRARY_PATH} REQUIRED)

  if(NOT GSL_INCLUDE_DIR OR NOT GSL_LIB OR NOT GSLCBLAS_LIB)
    message(FATAL_ERROR "Could not find GSL libraries or headers")
  else()
    message(STATUS "Found GSL include: ${GSL_INCLUDE_DIR}")
    message(STATUS "Found GSL lib: ${GSL_LIB}")
    message(STATUS "Found GSL CBLAS lib: ${GSLCBLAS_LIB}")
  endif()

  add_library(GSL::gsl UNKNOWN IMPORTED)
  set_target_properties(GSL::gsl PROPERTIES
    IMPORTED_LOCATION "${GSL_LIB}"
    INTERFACE_INCLUDE_DIRECTORIES "${GSL_INCLUDE_DIR}"
  )

  add_library(GSL::gslcblas UNKNOWN IMPORTED)
  set_target_properties(GSL::gslcblas PROPERTIES
    IMPORTED_LOCATION "${GSLCBLAS_LIB}"
    INTERFACE_INCLUDE_DIRECTORIES "${GSL_INCLUDE_DIR}"
  )

  # --------- Fetch Boost headers ---------
  set(BOOST_VERSION "boost-1.84.0")

  FetchContent_Declare(
    boost
    GIT_REPOSITORY https://github.com/boostorg/boost.git
    GIT_TAG ${BOOST_VERSION}
    GIT_SHALLOW TRUE
    GIT_SUBMODULES ""
  )

  FetchContent_MakeAvailable(boost)

  set(BOOST_ROOT ${boost_SOURCE_DIR})
  set(BOOST_MATH_INCLUDE_DIR ${BOOST_ROOT}/libs/math/include)
  set(BOOST_NUMERIC_ODEINT_INCLUDE_DIR ${BOOST_ROOT}/libs/numeric/odeint/include)

endif()

# ---- Add options ----
set(PYBHPT_MARCH native CACHE STRING "Value of the -march compiler option if supported by compiler")

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

# Include paths
if(DEFINED CMAKE_TOOLCHAIN_FILE)
  target_include_directories(cybhpt_full PRIVATE ${Boost_INCLUDE_DIRS})
else()
  target_include_directories(cybhpt_full PRIVATE
    ${BOOST_ROOT}
    ${BOOST_MATH_INCLUDE_DIR}
    ${BOOST_NUMERIC_ODEINT_INCLUDE_DIR}
  )
endif()

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
  target_compile_definitions(cybhpt_full PRIVATE PYBHPT_NO_OPENMP)
endif()

# Link libraries
target_link_libraries(cybhpt_full PRIVATE GSL::gsl GSL::gslcblas)

# Install (optional)
install(TARGETS cybhpt_full DESTINATION "./")
