cmake_minimum_required(VERSION 2.8)
project( PoKiTT CXX )

#-------------------------------------------------------------------------------
#-------------------------------- CONFIGURATION --------------------------------
#-------------------------------------------------------------------------------

option( ENABLE_TESTS "Enable tests" ON )

#-- Set the location for the installation.  Comment this out to install to /usr/local
if( CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT )
  set (CMAKE_INSTALL_PREFIX
    ${CMAKE_CURRENT_BINARY_DIR}/install/pokitt
    CACHE PATH "" FORCE
    )
endif( CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT )

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

find_package( Cantera )
set( TPL_LIBRARIES ${TPL_LIBRARIES} ${Cantera_LIBRARIES} )
set( TPL_INCLUDES ${TPL_INCLUDES} ${Cantera_INCLUDE_DIR} )


add_subdirectory( pokitt )

if( ENABLE_TESTS )
  if( NOT GIT_FOUND )
    message( SEND_ERROR "Git must be available to enable testing" )
  endif()
  
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