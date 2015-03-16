nebo_cuda_prep_dir()

file( GLOB XML_FILES ${PROJECT_SOURCE_DIR}/test/*.xml )
file( COPY ${XML_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

#-------------------------------------------------------------------------------
#--------------------------- Build examples ------------------------------------
#-------------------------------------------------------------------------------

include_directories( ${CMAKE_SOURCE_DIR} )
include_directories( ${CMAKE_BINARY_DIR}/config )

nebo_add_executable( rxndiff RxnDiffExample.cpp )
target_link_libraries( rxndiff pokitt )