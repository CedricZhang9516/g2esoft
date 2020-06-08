# CMakeLists.txt for TreeProc
# @author Taikan Suehara (Kyushu University)
cmake_minimum_required(VERSION 2.8)


#project( TreeProc ) 

aux_source_directory( ./src sources )
add_library( TreeProc SHARED ${sources} )
#target_include_directories( TreeProc PUBLIC "./include" )
target_link_libraries( TreeProc  ${ROOT_LIBRARIES} ${ROOT_COMPONENT_LIBRARIES} )

install( TARGETS TreeProc LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/lib )
