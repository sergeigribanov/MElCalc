cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(MElCalc)
find_package(MEUtils REQUIRED)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net Core Eve EG Geom Ged RGL)
include(${ROOT_USE_FILE})
file(GLOB HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h)
file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
add_library(SMEGen SHARED ${SOURCES})
target_include_directories(SMEGen
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE
  ${MEUtils_INCLUDE_DIRS})
target_link_libraries(SMEGen
  PRIVATE
  ${ROOT_LIBRARIES}
  ${MEUtils_LIBRARIES})
set(PUBLIC_HEADERS ${HEADERS})
set_target_properties(SMEGen
  PROPERTIES
  PUBLIC_HEADER "${PUBLIC_HEADERS}")
install(TARGETS SMEGen
		LIBRARY DESTINATION lib
		PUBLIC_HEADER DESTINATION include)

