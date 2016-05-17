install( FILES
  DiffusionCoeffMix.h
  SpeciesDiffusion.h
  ThermalCondMix.h
  ViscosityMix.h
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include/pokitt/transport
  PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
 )