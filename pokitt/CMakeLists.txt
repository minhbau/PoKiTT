add_subdirectory( kinetics )
add_subdirectory( thermo )
add_subdirectory( transport )

add_library( pokitt 
      CanteraObjects.cpp
      MixtureFraction.cpp
      kinetics/ReactionInfo.cpp
  )

target_link_libraries( pokitt
    PUBLIC
        exprlib
    PRIVATE
        ${Cantera_LIBRARIES}
        ${BLAS_LIBRARIES}
    )

target_include_directories( pokitt
    PUBLIC
        ${Cantera_INCLUDE_DIR}
        ${CMAKE_THREAD_LIBS_INIT}
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/config>
        $<INSTALL_INTERFACE:include>
    )

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
install( EXPORT PoKiTT
        FILE
            PoKiTT_Targets.cmake
        DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/pokitt
        )
