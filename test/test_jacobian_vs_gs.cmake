function( test_jacobian_vs_gs testName inputFile calcFile goldFile )

  add_test( NAME ${testName}_generate_output
            COMMAND aj_test --xml-input-file=${inputFile} --jacobian-vs-gs 
          )

  if( ${ARGC} EQUAL 4 )
    file( REMOVE_RECURSE ${CMAKE_CURRENT_BINARY_DIR}/${calcFile} )
    add_test( NAME ${testName}_compare
              COMMAND ${CMAKE_COMMAND} -E compare_files
              ${CMAKE_CURRENT_BINARY_DIR}/${calcFile}
              ${CMAKE_CURRENT_BINARY_DIR}/${goldFile}
            )
    set_property( TEST ${testName}_compare APPEND PROPERTY DEPENDS ${testName}_generate_output )
  endif( ${ARGC} EQUAL 4 )

endfunction()