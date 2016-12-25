cmake_minimum_required(VERSION 2.6)

Project(StructurePropagation)

# Required because of the c++0x used in Helpers
if(UNIX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
endif(UNIX)

set(EIGEN_DIR "none" CACHE FILEPATH "Path to Eigen 3.3.1 source")
include_directories(${EIGEN_DIR})

FIND_PACKAGE(VTK REQUIRED)
INCLUDE( ${VTK_USE_FILE} )

FIND_PACKAGE(ITK REQUIRED)
INCLUDE( ${ITK_USE_FILE} )

add_subdirectory(Mask)
add_subdirectory(DynamicProgramming)

ADD_EXECUTABLE(StructurePropagation StructurePropagation.h NonInteractive.cpp)
TARGET_LINK_LIBRARIES(StructurePropagation ${ITK_LIBRARIES} Mask)

