cmake_minimum_required(VERSION 3.1)
project( PoKiTT CXX C )

# warn about building in source
if( CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR )
  message( WARNING "\nIt appears that you are building directly in the source directory.\n"
    "This is strongly discouraged.\n"
    "Suggestion:\n\tmkdir build; cd build\n"
    "Then run cmake from the 'build' directory" ) 
endif( CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR )

set( CMAKE_CXX_STANDARD 11 )

if( ${CMAKE_HOST_APPLE} )
  message( "\nOn APPLE - adding -fvisibility-inlines-hidden \n")
  set( CMAKE_CXX_FLAGS "-fvisibility-inlines-hidden ${CMAKE_CXX_FLAGS}" )
endif()

#-------------------------------------------------------------------------------
#-------------------------------- CONFIGURATION --------------------------------
#-------------------------------------------------------------------------------

option( BUILD_UPSTREAM_LIBS "Auto-build ExprLib and SpatialOps (not recommended when linking to downstream apps)" OFF )
option( ENABLE_TESTS "Enable tests" ON )

#still a messy work in progress but it should compile and run
option( ENABLE_EXAMPLES "Enable reaction-diffusion and continuation method examples" ON )

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
  
  set( BOOST_INCLUDE CACHE PATH "Path to boost include directory" )
  set( NTHREADS 1 CACHE STRING "Number of threads if ENABLE_THREADS is ON" )

  # Build SpatialOps
  set( TPL_DIR ${PROJECT_BINARY_DIR}/tpl )
  file( MAKE_DIRECTORY ${TPL_DIR}/build/spatialops )
  execute_process( COMMAND ${GIT_EXECUTABLE} clone https://gitlab.multiscale.utah.edu/common/SpatialOps.git
                   WORKING_DIRECTORY ${TPL_DIR}
                   RESULT_VARIABLE result  )
  execute_process( COMMAND ${GIT_EXECUTABLE} reset --hard 117652f3bdbc2c9fb68042180772401a08ad06d2
                   WORKING_DIRECTORY ${TPL_DIR}/SpatialOps )

  execute_process( COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                   -DENABLE_THREADS=${ENABLE_THREADS}
                   -DENABLE_CUDA=${ENABLE_CUDA}
                   -DENABLE_TESTS=OFF
                   -DENABLE_EXAMPLES=OFF
                   -DCMAKE_INSTALL_PREFIX=${TPL_DIR}
                   -DBOOST_INCLUDEDIR=${BOOST_INCLUDE}
                   ${TPL_DIR}/SpatialOps
                   WORKING_DIRECTORY ${TPL_DIR}/build/spatialops
                   RESULT_VARIABLE result  )
  execute_process( COMMAND make -j4 install
                   WORKING_DIRECTORY ${TPL_DIR}/build/spatialops )
  
  # Build ExprLib
  file( MAKE_DIRECTORY ${TPL_DIR}/build/exprlib )
  execute_process( COMMAND ${GIT_EXECUTABLE} clone https://gitlab.multiscale.utah.edu/common/ExprLib.git
                   WORKING_DIRECTORY ${TPL_DIR}
                   RESULT_VARIABLE result )
  execute_process( COMMAND ${GIT_EXECUTABLE} reset --hard 89392f330a38892396a3578367294e4e2eaee2cb
                   WORKING_DIRECTORY ${TPL_DIR}/ExprLib )
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

find_library( Cantera_LIBRARIES NAMES cantera PATHS ${Cantera_DIR}/lib NO_DEFAULT_PATH )
find_path( Cantera_INCLUDE_DIR NAMES IdealGasMix.h PATHS ${Cantera_DIR}/include NO_DEFAULT_PATH PATH_SUFFIXES cantera )
if( ${Cantera_INCLUDE_DIR} STREQUAL Cantera_INCLUDE_DIR-NOTFOUND )
  message( WARNING "Could not find Cantera in suggested location ${Cantera_DIR}. Looking in other locations\n" )
  # try other locations
  find_library( Cantera_LIBRARIES NAMES cantera HINTS /usr/local )
  find_path( Cantera_INCLUDE_DIR NAMES IdealGasMix.h HINTS /usr/local PATH_SUFFIXES cantera )
  if( ${Cantera_INCLUDE_DIR} STREQUAL Cantera_INCLUDE_DIR-NOTFOUND )
    message( FATAL_ERROR "Could not find Cantera\n" )
  endif()
endif()
if( NOT ${Cantera_DIR} STREQUAL "" )
  string( REPLACE "include/cantera" "include" Cantera_INCLUDE_DIR ${Cantera_INCLUDE_DIR} )
endif()
message( STATUS "Cantera libraries at " ${Cantera_LIBRARIES} )
message( STATUS "Cantera includes  at " ${Cantera_INCLUDE_DIR} )


#-------------------------------------------------------------------------------
#------------------------------- BLAS & Threads --------------------------------
#-------------------------------------------------------------------------------
find_package( BLAS REQUIRED )

set( TPL_LIBRARIES ${TPL_LIBRARIES} ${Cantera_LIBRARIES}  ${BLAS_LIBRARIES} )
set( TPL_INCLUDES ${TPL_INCLUDES} ${Cantera_INCLUDE_DIR} )

# jcs it appears that on linux, the cantera library requires linkage with thread libraries.
find_package( Threads REQUIRED )
set( TPL_LIBRARIES ${TPL_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} )


#-------------------------------------------------------------------------------
#------------------------------ Build Everything -------------------------------
#-------------------------------------------------------------------------------

include_directories( 
  ${TPL_INCLUDES}
  ${CMAKE_SOURCE_DIR}
  )

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
