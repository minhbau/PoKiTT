
set( TPL_DIR ${PROJECT_BINARY_DIR}/tpl )

# Build SpatialOps
file( MAKE_DIRECTORY ${TPL_DIR}/build/spatialops )
execute_process( COMMAND ${GIT_EXECUTABLE} clone --depth 1 git://software.crsim.utah.edu/SpatialOps.git
                 WORKING_DIRECTORY ${TPL_DIR}
                 RESULT_VARIABLE result  )
execute_process( COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DENABLE_TESTS=OFF -DENABLE_THREADS=OFF -DENABLE_EXAMPLES=OFF -DCMAKE_INSTALL_PREFIX=${TPL_DIR} ${TPL_DIR}/SpatialOps
                 WORKING_DIRECTORY ${TPL_DIR}/build/spatialops
                 RESULT_VARIABLE result  )
execute_process( COMMAND make -j4 install
                 WORKING_DIRECTORY ${TPL_DIR}/build/spatialops )

# Build ExprLib
file( MAKE_DIRECTORY ${TPL_DIR}/build/exprlib )
execute_process( COMMAND ${GIT_EXECUTABLE} clone --depth 1 git://software.crsim.utah.edu/ExprLib.git
                 WORKING_DIRECTORY ${TPL_DIR}
                 RESULT_VARIABLE result )
execute_process( COMMAND ${CMAKE_COMMAND} -DENABLE_TESTS=OFF -DSpatialOps_DIR=${TPL_DIR}/share -DCMAKE_INSTALL_PREFIX=${TPL_DIR} -DENABLE_OUTPUT=OFF -DENABLE_LBMS=ON -DBUILD_GUI=OFF ${TPL_DIR}/ExprLib
                 WORKING_DIRECTORY ${TPL_DIR}/build/exprlib
                 RESULT_VARIABLE result )
execute_process( COMMAND make -j4 install
                 WORKING_DIRECTORY ${TPL_DIR}/build/exprlib )

find_package( ExprLib REQUIRED PATHS ${TPL_DIR}/share )


#-------------------------------------------------------------------------------
#--------------------------- Build Test executables ----------------------------
#-------------------------------------------------------------------------------

include_directories( ${CMAKE_SOURCE_DIR}/src )

#add_executable( diff_coef_test DiffusionCoeff_test.cpp )
#target_link_libraries( diff_coef_test exprlib )

#add_executable( enthalpy_test Enthalpy_test.cpp )
#target_link_libraries( enthalpy_test exprlib )

add_executable( cp_test HeatCapacity_Cp_test.cpp )
target_link_libraries( cp_test exprlib )

#add_executable( cv_test HeatCapacity_Cv_test.cpp )
#target_link_libraries( cv_test exprlib )

#add_executable( rxn_test ReactionRates_test.cpp )
#target_link_libraries( rxn_test exprlib )

#add_executable( temperature_test Temperature_test.cpp )
#target_link_libraries( temperature_test exprlib )

#add_executable( therm_cond_test ThermalCondMix_test.cpp )
#target_link_libraries( therm_cond_test exprlib )

#add_executable( viscosity_test ViscosityMix_test.cpp )
#target_link_libraries( viscosity_test exprlib )