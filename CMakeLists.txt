##
# our functions and macros
##
    # logging_level_to_define()
    # Gives the define corresponding to a global logging level.
    #
    # For example: The logging level 'INFO' corresponds to the define 'logINFO'.
    #
    # Parameters:
    #   LOGGING_LEVEL   STRING  One of the logging levels:
    #                           NO_LOGGING, ERROR, WARNING, INFO, DEBUG, DEBUG1, (...), DEBUG4
    #                           input parameter
    #
    #   DEFINE          STRING  output parameter
    # 
    MACRO(LOGGING_LEVEL_TO_DEFINE LOGGING_LEVEL DEFINE)
        IF(${LOGGING_LEVEL} STREQUAL "NO_LOGGING")
            SET(${DEFINE} "psf::logNO_LOGGING")
        ELSEIF(${LOGGING_LEVEL} STREQUAL "ERROR")
            SET(${DEFINE} "psf::logERROR")
        ELSEIF(${LOGGING_LEVEL} STREQUAL "WARNING")
            SET(${DEFINE} "psf::logWARNING")
        ELSEIF(${LOGGING_LEVEL} STREQUAL "INFO")
            SET(${DEFINE} "psf::logINFO")
        ELSEIF(${LOGGING_LEVEL} STREQUAL "DEBUG")
            SET(${DEFINE} "psf::logDEBUG")
        ELSEIF(${LOGGING_LEVEL} STREQUAL "DEBUG1")
            SET(${DEFINE} "psf::logDEBUG1")
        ELSEIF(${LOGGING_LEVEL} STREQUAL "DEBUG2")
            SET(${DEFINE} "psf::logDEBUG2")
        ELSEIF(${LOGGING_LEVEL} STREQUAL "DEBUG3")
            SET(${DEFINE} "psf::logDEBUG3")
        ELSEIF(${LOGGING_LEVEL} STREQUAL "DEBUG4")
            SET(${DEFINE} "psf::logDEBUG4")
        ELSE(${LOGGING_LEVEL} STREQUAL "NO_LOGGING")
            MESSAGE(SEND_ERROR "Unknown LOGGING_LEVEL: ${LOGGING_LEVEL}. Default to INFO.")
            SET(${DEFINE} "psf::logINFO")
        ENDIF(${LOGGING_LEVEL} STREQUAL "NO_LOGGING")
    ENDMACRO(LOGGING_LEVEL_TO_DEFINE)



##
# require CMake 2.6
##
    CMAKE_MINIMUM_REQUIRED(VERSION 2.6)


   
##
# prerequisites
##
    PROJECT(PSF)

    # setup testing
    ENABLE_TESTING()
    INCLUDE(CTest) 

    # toggle normal/advanced parameters
    MARK_AS_ADVANCED(CLEAR CMAKE_VERBOSE_MAKEFILE)

    IF(NOT CMAKE_BUILD_TYPE)
	    SET(CMAKE_BUILD_TYPE Release CACHE STRING
	        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
	        FORCE)
    ENDIF(NOT CMAKE_BUILD_TYPE)

    CONFIGURE_FILE(
        ${CMAKE_CURRENT_SOURCE_DIR}/cmake/config.h.cmake 
        ${PSF_BINARY_DIR}/include/psf/config.h
    )



##
# build options
##
    # logging level
    SET(LOGGING_LEVEL "INFO" CACHE STRING "Choose a global logging level: NO_LOGGING, ERROR, WARNING, INFO, DEBUG, DEBUG1, ..., DEBUG4")

    # vigra include
    SET(VIGRA_INCLUDE "/usr/include/" CACHE STRING "Path to vigra include dir.") 

    # build static or shared libraries
    OPTION(BUILD_SHARED_LIBS "Build psf and helper libraries in shared mode (else static)." ON)

### check host system type
# check for 64 bit OS
# Pointer has 8 bit on a 64Bit OS(only for intel&AMD)
    IF( NOT APPLE )
        IF( CMAKE_SIZEOF_VOID_P MATCHES 4 )
            SET(X86_64 FALSE)
        ELSE( CMAKE_SIZEOF_VOID_P MATCHES 4 )
            SET(X86_64 TRUE)
        ENDIF( CMAKE_SIZEOF_VOID_P MATCHES 4 )
    ENDIF( NOT APPLE )


##
# cflags and compiler switches
##
    IF(CMAKE_BUILD_TYPE STREQUAL "Debug")
        ADD_DEFINITIONS(-DPSFDEBUG)
    ENDIF(CMAKE_BUILD_TYPE STREQUAL "Debug")
    
    IF(MSVC)
	    ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE -D_SCL_SECURE_NO_WARNINGS -DEXP_STL)
    ELSE(MSVC)
	    ADD_DEFINITIONS(-Wall)
	    SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -D_FILE_OFFSET_BITS=64")
	    SET(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG -D_FILE_OFFSET_BITS=64")
	    SET(CMAKE_CXX_FLAGS_DEBUG  "-O0 -Werror -ggdb3 -D_FILE_OFFSET_BITS=64")
	    IF(WITH_GCOV AND CMAKE_BUILD_TYPE STREQUAL "Debug")
		    SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage")
	    ENDIF(WITH_GCOV AND CMAKE_BUILD_TYPE STREQUAL "Debug")
    ENDIF(MSVC)

##
# global logging level
#
    LOGGING_LEVEL_TO_DEFINE(LOGGING_LEVEL LOG_DEFINE)
    ADD_DEFINITIONS(-D FILELOG_MAX_LEVEL=${LOG_DEFINE})


##
# 64-bit support
##
    IF(X86_64)
	    ADD_DEFINITIONS(-DX86_64)
    ENDIF(X86_64)

##
# includes paths
##
    INCLUDE_DIRECTORIES(
	    ${PSF_BINARY_DIR}
        ${PSF_BINARY_DIR}/include
	    ${PSF_SOURCE_DIR}
	    ${PSF_SOURCE_DIR}/src
	    ${PSF_SOURCE_DIR}/include
        ${VIGRA_INCLUDE}
    )


##
# directories to build
##
    ADD_SUBDIRECTORY(src)

    IF(BUILD_TESTING) # user option
        ADD_SUBDIRECTORY(tests)
    ENDIF(BUILD_TESTING)

##
# summary
##
    MESSAGE(STATUS)
    MESSAGE(STATUS "-----------------------------------------")
    MESSAGE(STATUS "Running on: ${CMAKE_SYSTEM}")
    MESSAGE(STATUS "Running on: ${CMAKE_SYSTEM_VERSION}")
    MESSAGE(STATUS "This is a 64bit system: ${X86_64}")
    MESSAGE(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
    MESSAGE(STATUS "Global logging level: ${LOGGING_LEVEL}")
    MESSAGE(STATUS "-----------------------------------------")
    MESSAGE(STATUS)
