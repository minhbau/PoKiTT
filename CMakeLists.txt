cmake_minimum_required(VERSION 2.8)
project( PoKiTT CXX )

#-------------------------------------------------------------------------------
#-------------------------------- CONFIGURATION --------------------------------
#-------------------------------------------------------------------------------

option( BUILD_UPSTREAM_LIBS "Auto-build ExprLib and SpatialOps" OFF )
option( ENABLE_TESTS "Enable tests" ON )

#-- Set the location for the installation.  Comment this out to install to /usr/local
if( CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT )
  set (CMAKE_INSTALL_PREFIX
    ${CMAKE_CURRENT_BINARY_DIR}/install/pokitt
    CACHE PATH "" FORCE
    )
endif( CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT )

#-------------------------------------------------------------------------------
#------------------------------------- GIT -------------------------------------
#-------------------------------------------------------------------------------

# look for git.  This is used to configure version information into the
# executable and also to build upstream dependencies if necessary
find_package( Git )

if( GIT_FOUND )
  execute_process(
    COMMAND ${GIT_EXECUTABLE} log -1 "--pretty=format:\"%H\""
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE PoKiTT_REPO_HASH 
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  execute_process(
    COMMAND ${GIT_EXECUTABLE} log -1 "--pretty=format:\"%cd\""
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE PoKiTT_REPO_DATE 
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
else( GIT_FOUND )
  set( PoKiTT_REPO_DATE "\"No information available\"" )
  set( PoKiTT_REPO_HASH "\"No information available\"" )
endif( GIT_FOUND )

configure_file(
  PoKiTTVersion.h.in
  ${PROJECT_BINARY_DIR}/config/PoKiTTVersion.h
  @ONLY
  )


#-------------------------------------------------------------------------------
#----------------------------- ExprLib & SpatialOps ----------------------------
#-------------------------------------------------------------------------------

# NOTES:
#  If we build upstream libraries, we are building the "dumb" version: no
#  multithreading, no GPU.
#
#  For more general build configurations, you should build custom versions of
#  ExprLib and SpatialOps and point to those by setting:
#     ExprLib_DIR
#  to the location where ExprLib's configuration was installed (nominally install/share)
#
if( BUILD_UPSTREAM_LIBS )

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
  execute_process( COMMAND ${CMAKE_COMMAND} -DENABLE_TESTS=OFF -DSpatialOps_DIR=${TPL_DIR}/share -DCMAKE_INSTALL_PREFIX=${TPL_DIR} -DENABLE_OUTPUT=OFF -DBUILD_GUI=OFF ${TPL_DIR}/ExprLib
                   WORKING_DIRECTORY ${TPL_DIR}/build/exprlib
                   RESULT_VARIABLE result )
  execute_process( COMMAND make -j4 install
                   WORKING_DIRECTORY ${TPL_DIR}/build/exprlib )
  
  find_package( ExprLib REQUIRED PATHS ${TPL_DIR}/share )
  unset( TPL_DIR )

else( BUILD_UPSTREAM_LIBS )

  find_package( ExprLib REQUIRED )

endif( BUILD_UPSTREAM_LIBS )

set( TPL_INCLUDES ${TPL_INCLUDES} ${ExprLib_INCLUDE_DIR} )
set( TPL_LIBRARIES ${TPL_LIBRARIES} exprlib )


#-------------------------------------------------------------------------------
#----------------------------------- Cantera -----------------------------------
#-------------------------------------------------------------------------------

find_package( Cantera )
set( TPL_LIBRARIES ${TPL_LIBRARIES} ${Cantera_LIBRARIES} )
set( TPL_INCLUDES ${TPL_INCLUDES} ${Cantera_INCLUDE_DIR} )


#-------------------------------------------------------------------------------
#------------------------------ Build Everything -------------------------------
#-------------------------------------------------------------------------------

include_directories( ${TPL_INCLUDES} )

add_subdirectory( pokitt )

if( ENABLE_TESTS )
  enable_testing()
  add_subdirectory( test )
endif( ENABLE_TESTS )


#-------------------------------------------------------------------------------
#-------------------------------- INSTALLATION ---------------------------------
#-------------------------------------------------------------------------------
configure_file(
  ${PROJECT_SOURCE_DIR}/PoKiTTConfig.cmake.in
  ${CMAKE_BINARY_DIR}/config/PoKiTTConfig.cmake
  @ONLY
  )

install( FILES
  ${CMAKE_BINARY_DIR}/config/PoKiTTConfig.cmake
  DESTINATION ${CMAKE_INSTALL_PREFIX}/share
  )

install( FILES
  ${PROJECT_BINARY_DIR}/config/PoKiTTVersion.h
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include/pokitt
  )