cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(MElCalculator)
option(MEUTILS_DIR "The path to the MatrixElementUtils installation.")
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "~/install/SimpleMatrixElementGenerator")
endif()
if (NOT MEUTILS_DIR)
  if (DEFINED ENV{MEUTILS_DIR})
    set(MEUTILS_DIR $ENV{MEUTILS_DIR})
  else()
    set(MEUTILS_DIR "~/install/MatrixElementUtils")
  endif()
endif()
message("-- INSTALL PREFIX: " ${CMAKE_INSTALL_PREFIX})
message("-- The path to the MatrixElementUtils installation: " ${MEUTILS_DIR})
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net Core Eve EG Geom Ged RGL)
include(${ROOT_USE_FILE})
file(GLOB HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/inc/*.h)
file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
add_library(SMEGen SHARED ${SOURCES})
target_include_directories(SMEGen PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/inc PRIVATE ${MEUTILS_DIR}/inc)
target_link_libraries(SMEGen ${ROOT_LIBRARIES} ${MEUTILS_DIR}/lib/libMEUtils.so)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/sme-generator.sh.in sme-generator.sh)
set(PUBLIC_HEADERS ${HEADERS})
set_target_properties(SMEGen
  PROPERTIES
  PUBLIC_HEADER "${PUBLIC_HEADERS}"
  RESOURCE "${PROJECT_BINARY_DIR}/sme-generator.sh")
install(TARGETS SMEGen
		LIBRARY DESTINATION lib
		PUBLIC_HEADER DESTINATION inc
		RESOURCE DESTINATION env COMPONENT base)

