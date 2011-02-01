cmake_minimum_required(VERSION 2.6)

Project(StructurePropagation)

FIND_PACKAGE(Qt4 REQUIRED)
INCLUDE(${QT_USE_FILE})

QT4_WRAP_UI(UISrcs form.ui)
QT4_WRAP_CPP(MOCSrcs form.h ProgressThread.h)

include_directories(${include_directories} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

FIND_PACKAGE(VTK REQUIRED)
INCLUDE( ${USE_VTK_FILE} )

FIND_PACKAGE(ITK REQUIRED)
INCLUDE( ${USE_ITK_FILE} )

# OpenGM
SET(OpenGM_IncludePath "/media/portable/src/opengm10/include" CACHE STRING "Path to OpenGM include/ directory")
include_directories(${OpenGM_IncludePath})
FILE(GLOB OpenGM_Headers ${OpenGM_IncludePath}/opengm/* ${OpenGM_IncludePath}/opengm/inference/*)


ADD_EXECUTABLE(StructurePropagation main.cpp form.cpp ProgressThread.cpp
${OpenGM_Headers}
vtkScribbleInteractorStyle.cxx Helpers.cpp StructurePropagationBase.cpp
${MOCSrcs} ${UISrcs})
TARGET_LINK_LIBRARIES(StructurePropagation vtkHybrid vtkWidgets
QVTK ${QT_LIBRARIES}
ITKCommon ITKIO ITKStatistics
)
