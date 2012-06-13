/*
Copyright (C) 2011 David Doria, daviddoria@gmail.com

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

// ITK
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkVectorImage.h"

// Submodules
#include "Mask/Mask.h"
#include "Mask/ITKHelpers/ITKHelpers.h"

#include "StructurePropagation.h"

int main(int argc, char* argv[])
{
  if(argc != 4)
    {
    std::cerr << "Required: image mask patchRadius" << std::endl;
    exit(-1);
    }

  std::stringstream ss;
  for(unsigned int i = 1; i < argc; ++i)
  {
    ss << argv[i] << " ";
  }
  std::string imageFilename;
  std::string maskFilename;
  unsigned int patchRadius;

  ss >> imageFilename >> maskFilename >> patchRadius;

  std::cout << "Image: " << imageFilename << std::endl;
  std::cout << "Mask: " << maskFilename << std::endl;
  std::cout << "Patch radius: " << patchRadius << std::endl;

  typedef ITKHelpersTypes::UnsignedCharVectorImageType ImageType;
  ImageType::Pointer image = ImageType::New();
  ITKHelpers::ReadImage(imageFilename, image.GetPointer());

  Mask::Pointer mask = Mask::New();
  mask->Read(maskFilename);

  StructurePropagation<ImageType> structurePropagation;

  return EXIT_SUCCESS;
}

