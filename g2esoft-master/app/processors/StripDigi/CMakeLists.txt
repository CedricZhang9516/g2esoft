# CMakeLists.txt for TreeProc
# @author Taikan Suehara (Kyushu University)
cmake_minimum_required(VERSION 2.8)

get_filename_component(curdir ${CMAKE_CURRENT_SOURCE_DIR} NAME)

aux_source_directory( ./src sources )
add_library( ${curdir} SHARED ${sources} )
target_include_directories( ${curdir} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include" )
target_link_libraries( ${curdir} ${ROOT_LIBRARIES} ${ROOT_COMPONENT_LIBRARIES} )

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data)

file(GLOB DataFiles ${CMAKE_CURRENT_SOURCE_DIR}/data/*)
foreach(DataFile ${DataFiles})
  message("Copying data file ${DataFile} into ${CMAKE_CURRENT_BINARY_DIR}/data")
  configure_file(${DataFile} ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
endforeach()

install( TARGETS ${curdir} LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/lib )
