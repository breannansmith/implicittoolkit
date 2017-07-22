include( CMakeSourceFiles.txt )

add_executable( implicittoolkit_qt4 ${Headers} ${Sources} ${HWQ_MOC_CPP} )

set_target_properties( implicittoolkit_qt4 PROPERTIES AUTOMOC TRUE )

find_package( OpenGL REQUIRED )
target_include_directories( implicitutilities SYSTEM PRIVATE ${OPENGL_INCLUDE_DIR} )

find_package( Qt4 COMPONENTS QtCore QtGui QtOpenGL REQUIRED )
include( ${QT_USE_FILE} )

target_link_libraries( implicittoolkit_qt4 ${QT_LIBRARIES} ${OPENGL_LIBRARIES} implicitutilities marchingcubes )

execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/assets ${CMAKE_CURRENT_BINARY_DIR}/assets )
