nebo_cuda_prep_dir()

file( GLOB XML_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.xml )
file( COPY ${XML_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

#-------------------------------------------------------------------------------
#--------------------------- Build Test executables ----------------------------
#-------------------------------------------------------------------------------

include_directories( ${CMAKE_SOURCE_DIR} )

nebo_add_executable( idealgas_test IdealGas_test.cpp )
target_link_libraries( idealgas_test pokitt )

nebo_add_executable( rxn_test ReactionRates_test.cpp )
target_link_libraries( rxn_test pokitt )

nebo_add_executable( temperature_test Temperature_test.cpp )
target_link_libraries( temperature_test pokitt )

nebo_add_executable( thermo_test Thermo_test.cpp )
target_link_libraries( thermo_test pokitt )

nebo_add_executable( transport_test Transport_test.cpp )
target_link_libraries( transport_test pokitt )

add_executable( mixfrac_test mixfracTest.cpp )
target_link_libraries( mixfrac_test pokitt )

add_test( NAME thermo_test_const         COMMAND thermo_test       --xml-input-file=poly_tester.xml --phase=const         )
add_test( NAME thermo_test_nasa          COMMAND thermo_test       --xml-input-file=poly_tester.xml --phase=nasa          )
add_test( NAME thermo_test_shomate       COMMAND thermo_test       --xml-input-file=poly_tester.xml --phase=shomate       )

add_test( NAME thermo_test_mix_nas_sho   COMMAND thermo_test --mix --xml-input-file=poly_tester.xml --phase=nasa_shomate  )
add_test( NAME thermo_test_mix_nas_con   COMMAND thermo_test --mix --xml-input-file=poly_tester.xml --phase=nasa_const    )
add_test( NAME thermo_test_mix_sho_con   COMMAND thermo_test --mix --xml-input-file=poly_tester.xml --phase=shomate_const )

add_test( NAME temperature_test_nas_sho  COMMAND temperature_test  --xml-input-file=poly_tester.xml --phase=nasa_shomate  )
add_test( NAME temperature_test_nas_con  COMMAND temperature_test  --xml-input-file=poly_tester.xml --phase=nasa_const    )
add_test( NAME temperature_test_sho_con  COMMAND temperature_test  --xml-input-file=poly_tester.xml --phase=shomate_const )

add_test( NAME idealgas_test_h2          COMMAND idealgas_test     --xml-input-file=h2o2.xml     )
add_test( NAME rxn_test_h2               COMMAND rxn_test          --xml-input-file=h2o2.xml     )
add_test( NAME temperature_test_h2       COMMAND temperature_test  --xml-input-file=h2o2.xml     )
add_test( NAME thermo_test_mix_h2        COMMAND thermo_test --mix --xml-input-file=h2o2.xml     )
add_test( NAME thermo_test_h2            COMMAND thermo_test       --xml-input-file=h2o2.xml     )
add_test( NAME transport_test_h2         COMMAND transport_test    --xml-input-file=h2o2.xml     )

add_test( NAME idealgas_test_GRI         COMMAND idealgas_test     --xml-input-file=gri30.xml    )
add_test( NAME rxn_test_GRI              COMMAND rxn_test          --xml-input-file=gri30.xml    )
add_test( NAME temperature_test_GRI      COMMAND temperature_test  --xml-input-file=gri30.xml    )
add_test( NAME thermo_test_mix_GRI       COMMAND thermo_test --mix --xml-input-file=gri30.xml    )
add_test( NAME thermo_test_GRI           COMMAND thermo_test       --xml-input-file=gri30.xml    )
add_test( NAME transport_test_GRI        COMMAND transport_test    --xml-input-file=gri30.xml    )

add_test( NAME idealgas_test_methanol    COMMAND idealgas_test     --xml-input-file=methanol.xml )
add_test( NAME rxn_test_methanol         COMMAND rxn_test          --xml-input-file=methanol.xml )
add_test( NAME temperature_test_methanol COMMAND temperature_test  --xml-input-file=methanol.xml )
add_test( NAME thermo_test_mix_methanol  COMMAND thermo_test --mix --xml-input-file=methanol.xml )
add_test( NAME thermo_test_methanol      COMMAND thermo_test       --xml-input-file=methanol.xml )
add_test( NAME transport_test_methanol   COMMAND transport_test    --xml-input-file=methanol.xml )

add_test( NAME idealgas_test_ethanol     COMMAND idealgas_test     --xml-input-file=ethanol.xml  )
add_test( NAME rxn_test_ethanol          COMMAND rxn_test          --xml-input-file=ethanol.xml  )
add_test( NAME temperature_test_ethanol  COMMAND temperature_test  --xml-input-file=ethanol.xml  )
add_test( NAME thermo_test_mix_ethanol   COMMAND thermo_test --mix --xml-input-file=ethanol.xml  )
add_test( NAME thermo_test_ethanol       COMMAND thermo_test       --xml-input-file=ethanol.xml  )
add_test( NAME transport_test_ethanol    COMMAND transport_test    --xml-input-file=ethanol.xml  )

add_test( NAME mixfrac_test COMMAND mixfrac_test )