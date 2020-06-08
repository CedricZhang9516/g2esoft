# CMakeLists.txt for TreeProc
# @author Taikan Suehara (Kyushu University)
cmake_minimum_required(VERSION 2.8)

get_filename_component(curdir ${CMAKE_CURRENT_SOURCE_DIR} NAME)

aux_source_directory( ./src sources )
add_library( ${curdir} SHARED ${sources} )
target_include_directories( ${curdir} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include" )
target_link_libraries( ${curdir} ${ROOT_LIBRARIES} ${ROOT_COMPONENT_LIBRARIES} )

install( TARGETS ${curdir} LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/lib )
