About
-----
This software allows the user to specify a region to propagate.
This implementation is based on "Image Completion with Structure Propagation" by Jian Sun.

License
--------
GPLv3 (See LICENSE.txt)

Required dependencies
---------------------
- VTK 6.3
- ITK 4.10 (built with ITKV3_COMPATIBILITY=ON for itkImageToVectorImageFilter.h in ITKHelpers submodule (in ITKHelpers.hpp) )
- Qt 4.8.6
- Eigen 3.3.1

Build Instructions
------------------
- You must clone from github with:
git clone --recursive https://github.com/daviddoria/StructurePropagation

- Set QT_QMAKE_EXECUTABLE, ITK_DIR, VTK_DIR, and EIGEN_DIR (a custom CMake variable) via 'cmake -D' or using ccmake

Example
-------
From the build directory, run (where /path/to/data is the path to the 'data' directory in your clone of the StructurePropagation repository):
./StructurePropagation /path/to/data/trashcan.png /path/to/data/trashcan.mask /path/to/data/propagation.png 5 output.png

Notes
-----
- If you don't build with CMAKE_BUILD_TYPE=Release it will appear that the code is very very slow. With this set correctly, it should take ~5sec to run the above example on modern hardware
