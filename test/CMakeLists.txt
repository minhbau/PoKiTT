nebo_cuda_prep_dir()

file( GLOB XML_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.xml )
file( COPY ${XML_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

file( GLOB JACOBIAN_GOLD ${CMAKE_CURRENT_SOURCE_DIR}/jacobian-*_gold.dat )
file( COPY ${JACOBIAN_GOLD} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

#-------------------------------------------------------------------------------
#--------------------------- Build Test executables ----------------------------
#-------------------------------------------------------------------------------

include_directories( ${CMAKE_SOURCE_DIR} )

nebo_add_executable( idealgas_test IdealGas_test.cpp )
target_link_libraries( idealgas_test pokitt )

if( NOT ENABLE_CUDA )
    nebo_add_executable( idealgas_sens_test IdealGasSensitivity_test.cpp )
    target_link_libraries( idealgas_sens_test pokitt )
    
    nebo_add_executable( enthalpy_sens_test EnthalpySensitivity_test.cpp )
    target_link_libraries( enthalpy_sens_test pokitt )
endif()

nebo_add_executable( rxn_test ReactionRates_test.cpp )
target_link_libraries( rxn_test pokitt )

nebo_add_executable( aj_test AnalyticalJacobian_test.cpp )
target_link_libraries( aj_test pokitt )

nebo_add_executable( temperature_test Temperature_test.cpp )
target_link_libraries( temperature_test pokitt )

nebo_add_executable( thermo_test Thermo_test.cpp )
target_link_libraries( thermo_test pokitt )

nebo_add_executable( transport_test Transport_test.cpp )
target_link_libraries( transport_test pokitt )

add_executable( mixfrac_test mixfracTest.cpp )
target_link_libraries( mixfrac_test pokitt )

add_test( NAME thermo_test_const              COMMAND thermo_test            --xml-input-file=poly_tester.xml --phase=const         )
add_test( NAME thermo_test_nasa               COMMAND thermo_test            --xml-input-file=poly_tester.xml --phase=nasa          )
add_test( NAME thermo_test_mix_nas_con        COMMAND thermo_test --mix      --xml-input-file=poly_tester.xml --phase=nasa_const    )
                                                                             
add_test( NAME temperature_test_nas_con       COMMAND temperature_test       --xml-input-file=poly_tester.xml --phase=nasa_const    )
                                                                             
add_test( NAME idealgas_test_h2               COMMAND idealgas_test          --xml-input-file=h2o2.xml     )
add_test( NAME rxn_test_h2_1atm               COMMAND rxn_test               --xml-input-file=h2o2.xml --pressure=101325 )
add_test( NAME rxn_test_h2_2atm               COMMAND rxn_test               --xml-input-file=h2o2.xml --pressure=202650 )
add_test( NAME temperature_test_h2            COMMAND temperature_test       --xml-input-file=h2o2.xml     )
add_test( NAME thermo_test_mix_h2             COMMAND thermo_test --mix      --xml-input-file=h2o2.xml     )
add_test( NAME thermo_test_h2                 COMMAND thermo_test            --xml-input-file=h2o2.xml     )
add_test( NAME transport_test_h2              COMMAND transport_test         --xml-input-file=h2o2.xml     )
add_test( NAME ajrates_test_h2_1atm           COMMAND aj_test                --xml-input-file=h2o2.xml --pressure=101325 --rates-vs-cantera )
add_test( NAME ajrates_test_h2_2atm           COMMAND aj_test                --xml-input-file=h2o2.xml --pressure=202650 --rates-vs-cantera )  

if( NOT ENABLE_CUDA )
    add_test( NAME idealgas_sens_test_h2          COMMAND idealgas_sens_test )
    add_test( NAME species_enthalpy_sens_test_h2  COMMAND enthalpy_sens_test )
endif()
                                                                             
add_test( NAME idealgas_test_GRI              COMMAND idealgas_test          --xml-input-file=gri30.xml    )
add_test( NAME rxn_test_GRI_1atm              COMMAND rxn_test               --xml-input-file=gri30.xml --pressure=101325 )
add_test( NAME rxn_test_GRI_2atm              COMMAND rxn_test               --xml-input-file=gri30.xml --pressure=202650 )
add_test( NAME temperature_test_GRI           COMMAND temperature_test       --xml-input-file=gri30.xml    )
add_test( NAME thermo_test_mix_GRI            COMMAND thermo_test --mix      --xml-input-file=gri30.xml    )
add_test( NAME thermo_test_GRI                COMMAND thermo_test            --xml-input-file=gri30.xml    )
add_test( NAME transport_test_GRI             COMMAND transport_test         --xml-input-file=gri30.xml    )
if( CMAKE_BUILD_TYPE STREQUAL "Release" )
    add_test( NAME ajrates_test_GRI_1atm    COMMAND aj_test              --xml-input-file=gri30.xml --pressure=101325 --rates-vs-cantera )
    add_test( NAME ajrates_test_GRI_2atm    COMMAND aj_test              --xml-input-file=gri30.xml --pressure=202650 --rates-vs-cantera )
endif()
                                           
add_test( NAME idealgas_test_methanol         COMMAND idealgas_test          --xml-input-file=methanol.xml )
add_test( NAME rxn_test_methanol_1atm         COMMAND rxn_test               --xml-input-file=methanol.xml --pressure=101325 )
add_test( NAME rxn_test_methanol_2atm         COMMAND rxn_test               --xml-input-file=methanol.xml --pressure=202650 )
add_test( NAME temperature_test_methanol      COMMAND temperature_test       --xml-input-file=methanol.xml )
add_test( NAME thermo_test_mix_methanol       COMMAND thermo_test --mix      --xml-input-file=methanol.xml )
add_test( NAME thermo_test_methanol           COMMAND thermo_test            --xml-input-file=methanol.xml )
add_test( NAME transport_test_methanol        COMMAND transport_test         --xml-input-file=methanol.xml )
                                           
add_test( NAME idealgas_test_ethanol          COMMAND idealgas_test          --xml-input-file=ethanol.xml  )
add_test( NAME rxn_test_ethanol_1atm          COMMAND rxn_test               --xml-input-file=ethanol.xml --pressure=101325 )
add_test( NAME rxn_test_ethanol_2atm          COMMAND rxn_test               --xml-input-file=ethanol.xml --pressure=202650 )
add_test( NAME temperature_test_ethanol       COMMAND temperature_test       --xml-input-file=ethanol.xml  )
add_test( NAME thermo_test_mix_ethanol        COMMAND thermo_test --mix      --xml-input-file=ethanol.xml  )
add_test( NAME thermo_test_ethanol            COMMAND thermo_test            --xml-input-file=ethanol.xml  )
add_test( NAME transport_test_ethanol         COMMAND transport_test         --xml-input-file=ethanol.xml  )

add_test( NAME mixfrac_test COMMAND mixfrac_test )

include( test_jacobian_vs_gs.cmake )
test_jacobian_vs_gs( aj_h2o2      h2o2.xml      jacobian-h2o2.dat jacobian-h2o2_gold.dat )

if( CMAKE_BUILD_TYPE STREQUAL "Release" )
    test_jacobian_vs_gs( aj_gri30    gri30.xml    jacobian-gri30.dat    jacobian-gri30_gold.dat )
endif()








