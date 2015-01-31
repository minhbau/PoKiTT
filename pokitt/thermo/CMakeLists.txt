install( FILES
  Density.h
  Enthalpy.h
  InternalEnergy.h
  HeatCapacity_Cp.h
  HeatCapacity_Cv.h
  Pressure.h
  Temperature.h
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include/pokitt/thermo
  PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
 )