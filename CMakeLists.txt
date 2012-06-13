cmake_minimum_required(VERSION 2.6)

Project(StructurePropagation)

# Required because of the c++0x used in Helpers
if(UNIX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
endif(UNIX)

include_directories(/media/portable/src/Eigen)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE( ${USE_VTK_FILE} )

FIND_PACKAGE(ITK REQUIRED)
INCLUDE( ${USE_ITK_FILE} )

add_subdirectory(Mask)
add_subdirectory(DynamicProgramming)

ADD_EXECUTABLE(StructurePropagation NonInteractive.cpp)
TARGET_LINK_LIBRARIES(StructurePropagation ${ITK_LIBRARIES} Mask)

