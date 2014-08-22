nebo_cuda_prep_dir()

file( GLOB XML_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.xml )
file( COPY ${XML_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

#-------------------------------------------------------------------------------
#--------------------------- Build Test executables ----------------------------
#-------------------------------------------------------------------------------

include_directories( ${CMAKE_SOURCE_DIR} )

#nebo_add_executable( diff_coef_test DiffusionCoeff_test.cpp )
#target_link_libraries( diff_coef_test pokitt )
#add_test( diff_coef_test diff_coef_test )

nebo_add_executable( enthalpy_test Enthalpy_test.cpp )
target_link_libraries( enthalpy_test pokitt )
add_test( enthalpy_test enthalpy_test )

nebo_add_executable( cp_test HeatCapacity_Cp_test.cpp )
target_link_libraries( cp_test pokitt )
add_test( cp_test cp_test )

nebo_add_executable( cv_test HeatCapacity_Cv_test.cpp )
target_link_libraries( cv_test pokitt )
add_test( cv_test cv_test )

nebo_add_executable( rxn_test ReactionRates_test.cpp )
target_link_libraries( rxn_test pokitt )
add_test( rxn_test rxn_test )

nebo_add_executable( temperature_test Temperature_test.cpp )
target_link_libraries( temperature_test pokitt )
add_test( temperature_test temperature_test )

nebo_add_executable( therm_cond_test ThermalCondMix_test.cpp )
target_link_libraries( therm_cond_test pokitt )
add_test( therm_cond_test therm_cond_test )

nebo_add_executable( viscosity_test ViscosityMix_test.cpp )
target_link_libraries( viscosity_test pokitt )
add_test( viscosity_test viscosity_test )
