# CMakeLists.txt for TreeProc
# @author Taikan Suehara (Kyushu University)
cmake_minimum_required(VERSION 2.8)

set(libname g2esoftCommon)
set(dictname g2esoftDict)

# processors in subdirectory
add_subdirectory( processors )

# dictionary files
file(GLOB dict_include_files include/*.h)

# dictionary
include(${ROOT_USE_FILE})
ROOT_GENERATE_DICTIONARY(${dictname} ${dict_include_files} LINKDEF LinkDef.h)

aux_source_directory( ./src sources )
add_library( ${libname} SHARED ${sources} ${dictname}.cxx)
target_link_libraries( ${libname} ${ROOT_LIBRARIES} ${ROOT_COMPONENT_LIBRARIES} )

install( TARGETS ${libname} LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/lib )
install( FILES ${CMAKE_CURRENT_BINARY_DIR}/lib${dictname}_rdict.pcm DESTINATION ${CMAKE_SOURCE_DIR}/lib )
