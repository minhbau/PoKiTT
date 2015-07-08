cmake_minimum_required(VERSION 2.8)
project( PoKiTT CXX C )

#-------------------------------------------------------------------------------
#-------------------------------- CONFIGURATION --------------------------------
#-------------------------------------------------------------------------------

option( BUILD_UPSTREAM_LIBS "Auto-build ExprLib and SpatialOps (not recommended when linking to downstream apps)" OFF )
option( ENABLE_TESTS "Enable tests" ON )

#still a messy work in progress but it should compile and run
option( ENABLE_EXAMPLES "Enable reaction-diffusion example" ON )

set( Cantera_DIR "" CACHE PATH "Location where cantera was installed" )


# default to release builds
if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE "Release" CACHE STRING
       "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
       FORCE
       )
endif( NOT CMAKE_BUILD_TYPE )

#-- Set the location for the installation.  Comment this out to install to /usr/local
if( CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT )
  set( CMAKE_INSTALL_PREFIX
       ${CMAKE_CURRENT_BINARY_DIR}/install
       CACHE PATH "Installation directory" FORCE
       )
endif( CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT )

set( PoKiTT_CONFIG_INSTALL
     ${CMAKE_INSTALL_PREFIX}/share
     CACHE PATH
     "Location to install configuration files" )

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
    RESULT_VARIABLE RESULT
    )
    if( NOT ${RESULT} EQUAL 0 )
      set( PoKiTT_REPO_HASH "\"No information available\"" )
    endif()
  execute_process(
    COMMAND ${GIT_EXECUTABLE} log -1 "--pretty=format:\"%cd\""
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE PoKiTT_REPO_DATE 
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE RESULT
    )
    if( NOT ${RESULT} EQUAL 0 )
      set( PoKiTT_REPO_DATE "\"No information available\"" )
    endif()
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

  option( ENABLE_CUDA    "Enable GPU backend in SpatialOps"       OFF )
  option( ENABLE_THREADS "Enable multicore backend in SpatialOps" OFF )
  
  set( NTHREADS 1 CACHE STRING "Number of threads if ENABLE_THREADS is ON" )

  # Build SpatialOps
  set( TPL_DIR ${PROJECT_BINARY_DIR}/tpl )
  file( MAKE_DIRECTORY ${TPL_DIR}/build/spatialops )
  execute_process( COMMAND ${GIT_EXECUTABLE} clone --depth 1 https://software.crsim.utah.edu:8443/James_Research_Group/SpatialOps.git
                   WORKING_DIRECTORY ${TPL_DIR}
                   RESULT_VARIABLE result  )
  execute_process( COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                   -DENABLE_THREADS=${ENABLE_THREADS}
                   -DENABLE_CUDA=${ENABLE_CUDA}
                   -DENABLE_TESTS=OFF
                   -DENABLE_EXAMPLES=OFF
                   -DCMAKE_INSTALL_PREFIX=${TPL_DIR} 
                   ${TPL_DIR}/SpatialOps
                   WORKING_DIRECTORY ${TPL_DIR}/build/spatialops
                   RESULT_VARIABLE result  )
  execute_process( COMMAND make -j4 install
                   WORKING_DIRECTORY ${TPL_DIR}/build/spatialops )
  
  # Build ExprLib
  file( MAKE_DIRECTORY ${TPL_DIR}/build/exprlib )
  execute_process( COMMAND ${GIT_EXECUTABLE} clone --depth 1 https://software.crsim.utah.edu:8443/James_Research_Group/ExprLib.git
                   WORKING_DIRECTORY ${TPL_DIR}
                   RESULT_VARIABLE result )
  execute_process( COMMAND ${CMAKE_COMMAND}
                   -DENABLE_TESTS=OFF
                   -DSpatialOps_DIR=${TPL_DIR}/share
                   -DCMAKE_INSTALL_PREFIX=${TPL_DIR}
                   -DENABLE_OUTPUT=OFF
                   -DBUILD_GUI=OFF
                   ${TPL_DIR}/ExprLib
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
set( TPL_LIBRARIES ${TPL_LIBRARIES} ${Cantera_LIBRARIES}   )
set( TPL_INCLUDES  ${TPL_INCLUDES}  ${Cantera_INCLUDE_DIR} )


#-------------------------------------------------------------------------------
#------------------------------ Build Everything -------------------------------
#-------------------------------------------------------------------------------

include_directories( ${TPL_INCLUDES} )

add_subdirectory( pokitt )

if( ENABLE_TESTS )
  enable_testing()
  add_subdirectory( test )
endif( ENABLE_TESTS )

if( ENABLE_EXAMPLES )
  add_subdirectory( examples/ReactionDiffusion )
endif( ENABLE_EXAMPLES )


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
  DESTINATION ${PoKiTT_CONFIG_INSTALL}
  )

install( FILES
  ${PROJECT_BINARY_DIR}/config/PoKiTTVersion.h
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include/pokitt
  )