add_subdirectory( kinetics )
add_subdirectory( thermo )
add_subdirectory( transport )

add_library( pokitt CanteraObjects.cpp MixtureFraction.cpp )
target_link_libraries( pokitt ${TPL_LIBRARIES} )

install( FILES
  CanteraObjects.h
  MixtureFraction.h
  MixtureFractionExpr.h
  MixtureMolWeight.h
  MassToMoleFracs.h
  MoleToMassFracs.h
  SpecificToVolumetric.h
  VolumetricToSpecific.h
  SpecFromMixfrac.h
  SpeciesN.h
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include/pokitt
  PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
 )
 
install( TARGETS pokitt 
         EXPORT PoKiTT
         DESTINATION lib/pokitt )
install( EXPORT PoKiTT DESTINATION lib/pokitt )
