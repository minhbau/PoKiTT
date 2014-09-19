nebo_cuda_prep_dir()

file( GLOB XML_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.xml )
file( COPY ${XML_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

#-------------------------------------------------------------------------------
#--------------------------- Build Test executables ----------------------------
#-------------------------------------------------------------------------------

include_directories( ${CMAKE_SOURCE_DIR} )

nebo_add_executable( rxn_test ReactionRates_test.cpp )
target_link_libraries( rxn_test pokitt )
add_test( NAME rxn_test_h2      COMMAND rxn_test --xml-input-file=h2o2.xml         )
add_test( NAME rxn_test_GRI     COMMAND rxn_test --xml-input-file=gri30.xml        )
add_test( NAME rxn_test_ethanol COMMAND rxn_test --xml-input-file=ethanol_mech.xml )

nebo_add_executable( temperature_test Temperature_test.cpp )
target_link_libraries( temperature_test pokitt )
add_test( NAME temperature_test         COMMAND temperature_test --xml-input-file=poly_tester.xml)
add_test( NAME temperature_test_shomate COMMAND temperature_test --xml-input-file=poly_tester.xml --phase=shomate)
add_test( NAME temperature_test_GRI     COMMAND temperature_test --xml-input-file=gri30.xml)

nebo_add_executable( thermo_test Thermo_test.cpp )
target_link_libraries( thermo_test pokitt )
add_test( NAME thermo_test_mix          COMMAND thermo_test --xml-input-file=poly_tester.xml --mix )
add_test( NAME thermo_test              COMMAND thermo_test --xml-input-file=poly_tester.xml       )
add_test( NAME thermo_test_shomate_mix  COMMAND thermo_test --xml-input-file=poly_tester.xml --mix --phase=shomate )
add_test( NAME thermo_test_shomate      COMMAND thermo_test --xml-input-file=poly_tester.xml       --phase=shomate )
add_test( NAME thermo_test_GRI_mix      COMMAND thermo_test --xml-input-file=gri30.xml       --mix --phase=gri30   )

nebo_add_executable( transport_test Transport_test.cpp )
target_link_libraries( transport_test pokitt )
add_test( NAME transport_test COMMAND transport_test --xml-input-file=h2o2.xml )