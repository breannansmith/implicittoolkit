# CMake 3.6.1 fixes a bug with GCC and the system flag on OS X. Older CMake versions
# are acceptable if a different platform or toolchain is used.
cmake_minimum_required( VERSION 3.6.1 )

# Check that the Eigen is installed
execute_process( COMMAND ./scripts_include/eigen_found.sh RESULT_VARIABLE EIGEN_NOT_FOUND WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
if( EIGEN_NOT_FOUND )
  message( FATAL_ERROR "Failed to locate Eigen or checksum is incorrect. Please run get_dependencies.sh." )
endif()

project( ImplicitToolkit C CXX )

# If the user provided a CMAKE_BUILD_TYPE, create a cache variable
if( NOT CMAKE_BUILD_TYPE STREQUAL "" )
  set( CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel Coverage." FORCE )
# Otherwise default to Release
else( NOT CMAKE_BUILD_TYPE STREQUAL "" )
  set( CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel Coverage." FORCE )
endif( NOT CMAKE_BUILD_TYPE STREQUAL "" )

# Enforce the use of a valid build type
set( VALID_BUILD_TYPES None Debug Release RelWithDebInfo MinSizeRel Coverage )
if( NOT ";${VALID_BUILD_TYPES};" MATCHES ";${CMAKE_BUILD_TYPE};" )
  message( FATAL_ERROR "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel Coverage." )
endif()

# Disable compiler language extensions
set( CMAKE_CXX_EXTENSIONS OFF )
# Enable C++14
set( CMAKE_CXX_STANDARD 14 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )

# Append an extra space if flags are already set
if( CMAKE_CXX_FLAGS )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} " )
endif( CMAKE_CXX_FLAGS )

# Allow the user to optionally enable compiler sanitizers
set( SANITIZER "none" CACHE STRING "Enable compiler sanitizer modes, options are: none, address, undefined" )

# Enable compiler specific flags
if( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" )
  include( CMakeListsGNUCXXFlags.txt )
elseif( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel" )
  include( CMakeListsIntelCXXFlags.txt )
elseif( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" )
  include( CMakeListsClangCXXFlags.txt )
elseif( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang" )
  # For now, same options to Clang and AppleClang. May change in the future.
  include( CMakeListsClangCXXFlags.txt )
else()
  message( FATAL_ERROR "Unsupported C++ compiler ${CMAKE_CXX_COMPILER_ID}. Please email the maintainer." )
endif()

set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake )

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

option( STRICT_BUILD "Enable aggressive warnings and treat warnings as errors" OFF )
option( USE_QT4 "Build Qt4 front ends" ON )
option( USE_HDF5 "Serialization support with HDF5" ON )
option( USE_QHULL "Support for Qhull in the implicit toolkit" ON )

add_subdirectory( marchingcubes )
add_subdirectory( implicitutilities )
add_subdirectory( implicittoolkitcli )
if( USE_QT4 )
  add_subdirectory( implicittoolkitqt4 )
endif()
