include( CMakeSourceFiles.txt )

add_library( implicitutilities ${Headers} ${Sources} )

find_package( Eigen3 3.2.5 REQUIRED )
target_include_directories( implicitutilities SYSTEM PUBLIC ${EIGEN3_INCLUDE_DIR} )

# TODO: Fix the code so it can work with dynamically linked version of qhull
if( USE_QHULL )
  find_package( Qhull REQUIRED )
  if( NOT ${QHULL_FOUND} )
    message( FATAL_ERROR "QHull not found." )
  endif()
  message( STATUS "QHULL Directories: ${QHULL_INCLUDE_DIRS}" )
  target_include_directories( implicitutilities SYSTEM PRIVATE ${QHULL_INCLUDE_DIRS} )
  set_source_files_properties( Qhull.cpp PROPERTIES COMPILE_DEFINITIONS USE_QHULL )
  target_link_libraries( implicitutilities INTERFACE ${QHULL_LIBRARIES} )
endif()

if( USE_HDF5 )
  find_package( HDF5 REQUIRED )
  if( NOT ${HDF5_FOUND} )
    message( FATAL_ERROR "HDF5 not found." )
  endif()
  target_include_directories( implicitutilities SYSTEM PUBLIC ${HDF5_INCLUDE_DIRS} )
  target_compile_definitions( implicitutilities PUBLIC USE_HDF5 )
  target_link_libraries( implicitutilities INTERFACE ${HDF5_LIBRARIES} )
endif()
