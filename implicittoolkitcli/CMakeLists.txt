add_executable( implicittoolkit_cli "implicittoolkit_cli.cpp" "parse_params.h" )

#set_target_properties( implicittoolkit_cli PROPERTIES AUTOMOC TRUE )

target_link_libraries( implicittoolkit_cli implicitutilities marchingcubes )
