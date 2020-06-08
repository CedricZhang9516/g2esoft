# CMakeLists.txt for g2esoft
# @author Taikan Suehara (Kyushu University)
cmake_minimum_required(VERSION 2.8)

project( g2esoft ) 
set(exec_name ${PROJECT_NAME}.bin)

# root inclusion
find_package(ROOT 6.0 REQUIRED COMPONENTS XMLIO XMLParser)

# global include directories
include_directories( "include" )
include_directories( ${ROOT_INCLUDE_DIRS} )

# global compiler definitions
add_definitions( "-std=c++11" )
add_definitions( ${ROOT_DEFINITIONS} )

# sub projects
add_subdirectory( TreeProc )
add_subdirectory( app )

add_executable( ${exec_name} g2esoft.cxx )
target_link_libraries( ${exec_name} TreeProc ${ROOT_LIBRARIES} ${ROOT_COMPONENT_LIBRARIES} )

install( TARGETS ${exec_name} DESTINATION ${CMAKE_SOURCE_DIR}/bin )
