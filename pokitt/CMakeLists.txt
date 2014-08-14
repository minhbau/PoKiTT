add_subdirectory( kinetics )
add_subdirectory( thermo )
add_subdirectory( transport )

install( FILES
  CanteraObjects.h
  MixtureMoleWeightExpr.h
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include/pokitt
  PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
 )
 