# FindCantera.cmake
#
# Used to locate a cantera installation.
#
#  find_package( Cantera )
#
# Sets the following:
#  Cantera_FOUND = [true|false]
#  cantera - a target to link against

include(FindPackageHandleStandardArgs)
include(CMakePrintHelpers)

# some system tools to help locate the library
find_package( PkgConfig QUIET )
pkg_check_modules( PC_cantera
        QUIET
        cantera>=2.6
        cantera
    )

#cmake_print_variables(
#        PC_cantera_FOUND
#        PC_cantera_LIBRARIES
#        PC_cantera_LINK_LIBRARIES
#        PC_cantera_LIBRARY_DIRS
#        PC_cantera_INCLUDE_DIRS
#        PC_cantera_VERSION
#    )

find_path( Cantera_INCLUDE_DIR
        NAMES thermo.h
        PATHS ${Cantera_DIR} ${PC_cantera_INCLUDE_DIRS}
        PATH_SUFFIXES include include/cantera
        DOC "Path to Cantera header files"
    )

find_library( Cantera_LIBRARY
        NAMES cantera
        HINTS ${Cantera_DIR} ${PC_cantera_LIBRARY_DIRS}
        PATH_SUFFIXES lib
        NO_DEFAULT_PATH
    )

if( NOT ${Cantera_INCLUDE_DIR} STREQUAL Cantera_INCLUDE_DIR-NOTFOUND )

    string( REPLACE "include/cantera" "include" Cantera_INCLUDE_DIR ${Cantera_INCLUDE_DIR} )

    set( Cantera_config_file ${Cantera_INCLUDE_DIR}/cantera/base/config.h )

    ## Version information
    file( STRINGS ${Cantera_config_file}
            version_string
            LIMIT_COUNT 1000
            REGEX "(CANTERA_SHORT_VERSION)"
        )
    string( REGEX MATCH "([0-9]\.[0-9]+)" Cantera_VERSION  ${version_string} )
    string( REGEX MATCHALL "([0-9]+)" version_split ${Cantera_VERSION} )
    # Parse into separate numbers
    list( GET version_split 1 Cantera_VERSION_MINOR )
    list( GET version_split 0 Cantera_VERSION_MAJOR )

    # Check all required information
    find_package_handle_standard_args( Cantera
        FOUND_VAR Cantera_FOUND
        REQUIRED_VARS
            Cantera_INCLUDE_DIR
            Cantera_LIBRARY
        VERSION_VAR
            Cantera_VERSION
    )

    if( NOT TARGET cantera )
        add_library( cantera UNKNOWN IMPORTED )
        set_target_properties( cantera PROPERTIES
                IMPORTED_LOCATION "${Cantera_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${Cantera_INCLUDE_DIR}"
            )
    endif()

#    cmake_print_variables( Cantera_INCLUDE_DIR Cantera_LIBRARY Cantera_PREFIX Cantera_VERSION )
    mark_as_advanced(
        Cantera_INCLUDE_DIR
        Cantera_LIBRARY
        Cantera_VERSION
    )

    # See if we need to grab YAML as well
    file( STRINGS ${Cantera_config_file}
            external_yaml
            LIMIT_COUNT 1000
            REGEX "(#define CT_USE_SYSTEM_YAMLCPP)"
        )

    if( NOT external_yaml STREQUAL "" )

        pkg_check_modules( PC_yamlcpp
                QUIET
                yaml-cpp
            )
#        cmake_print_variables(
#                PC_yamlcpp_FOUND
#                PC_yamlcpp_LIBRARIES
#                PC_yamlcpp_LINK_LIBRARIES
#                PC_yamlcpp_LIBRARY_DIRS
#                PC_yamlcpp_INCLUDE_DIRS
#                PC_yamlcpp_PREFIX
#                PC_yamlcpp_VERSION
#            )

        find_package( yaml-cpp REQUIRED
                PATHS
                    ${PC_cantera_INCLUDE_DIRS}
                    ${YAML_DIR}
                PATH_SUFFIXES cmake/yaml-cpp
            )
#        cmake_print_properties( TARGETS yaml-cpp PROPERTIES LOCATION INTERFACE_INCLUDE_DIRECTORIES )
#        cmake_print_variables( yaml-cpp_LIBRARIES yaml-cpp_INCLUDE_DIR yaml-cpp_VERSION )

        find_package( BLAS REQUIRED )
        find_package( Threads REQUIRED )

        target_link_libraries( cantera INTERFACE yaml-cpp ${BLAS_LIBRARIES} Threads::Threads )

    endif()

endif()
