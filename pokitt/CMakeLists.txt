add_subdirectory( kinetics )
add_subdirectory( thermo )
add_subdirectory( transport )

add_library( pokitt CanteraObjects.cpp )
target_link_libraries( pokitt ${TPL_LIBRARIES} )

install( FILES
  CanteraObjects.h
  MixtureMolWeight.h
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include/pokitt
  PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
 )
 
install( TARGETS pokitt EXPORT pokitt DESTINATION lib/pokitt ) 