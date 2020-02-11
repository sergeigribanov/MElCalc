cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(MElCalc)
find_package(MEUtils REQUIRED)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net Core Eve EG Geom Ged RGL)
include(${ROOT_USE_FILE})
file(GLOB HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp)
file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
add_library(MElCalc SHARED ${SOURCES})
target_include_directories(MElCalc
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/MElCalc>
  PRIVATE
  ${MEUtils_INCLUDE_DIRS})
target_link_libraries(MElCalc
  PRIVATE
  ${ROOT_LIBRARIES}
  ${MEUtils_LIBRARIES})
set(PUBLIC_HEADERS ${HEADERS})
set_target_properties(MElCalc
  PROPERTIES
  PUBLIC_HEADER "${PUBLIC_HEADERS}")
install(TARGETS MElCalc
  EXPORT MElCalcTargets
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include/MElCalc)
install(EXPORT MElCalcTargets
  FILE MElCalcTargets.cmake
  DESTINATION lib/cmake)
include(CMakePackageConfigHelpers)
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/MElCalcConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/MElCalcConfig.cmake"
  INSTALL_DESTINATION "lib/cmake/example"
  NO_SET_AND_CHECK_MACRO
  NO_CHECK_REQUIRED_COMPONENTS_MACRO)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/MElCalcConfigVersion.cmake"
  VERSION "${MElCalc_VERSION_MAJOR}.${MElCalc_VERSION_MINOR}"
  COMPATIBILITY AnyNewerVersion)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/MElCalcConfig.cmake
  DESTINATION lib/cmake)
export(EXPORT MElCalcTargets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/MElCalcTargets.cmake")
export(PACKAGE MElCalc)
