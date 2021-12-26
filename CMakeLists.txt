cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

project(HairSimu VERSION 1.0)
message(STATUS "[Project] ${PROJECT_NAME}")

# Configuration types
set(CONFIGURATION_TYPES "Debug;Release;RelWithDebInfo")
if(win32)
  if(NOT OIDN_DEFAULT_CMAKE_CONFIGURATION_TYPES_SET)
    set(CMAKE_CONFIGURATION_TYPES "${CONFIGURATION_TYPES}"
        CACHE STRING "List of generated configurations." FORCE)
  endif()
else()
  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the build type." FORCE)
  endif()
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${CONFIGURATION_TYPES})
endif()

if(CMAKE_INSTALL_PREFIX_INITIALIXED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR} CACHE PATH "install prefix" FORCE)
endif()

set(CMAKE_DEBUG_POSTFIX d)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

file( COPY ${PROJECT_SOURCE_DIR}/world DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

set(FREEGLUT_DIR "./ext/freeglut-3.0.0")
INCLUDE_DIRECTORIES("${FREEGLUT_DIR}/include")
 
if(CMAKE_BUILD_TYPE MATCHES "Debug")
    link_directories("${FREEGLUT_DIR}/lib/Debug")
elseif(CMAKE_BUILD_TYPE MATCHES "Release")
    link_directories("${FREEGLUT_DIR}/lib/Release")
endif()

add_subdirectory(freeglut)
add_subdirectory(src)