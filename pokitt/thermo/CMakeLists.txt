install( FILES
  Enthalpy.h
  HeatCapacity_Cp.h
  HeatCapacity_Cv.h
  Temperature.h
  TemperaturePowers.h
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include/pokitt/thermo
  PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
 )