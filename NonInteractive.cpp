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
  if(argc != 6)
    {
    std::cerr << "Required: image mask propagationImage patchRadius output" << std::endl;
    exit(-1);
    }

  std::stringstream ss;
  for(unsigned int i = 1; i < argc; ++i)
  {
    ss << argv[i] << " ";
  }
  std::string imageFilename;
  std::string maskFilename;
  std::string propagationFilename;
  unsigned int patchRadius;
  std::string outputFilename;

  ss >> imageFilename >> maskFilename >> propagationFilename >> patchRadius >> outputFilename;

  std::cout << "Image: " << imageFilename << std::endl;
  std::cout << "Mask: " << maskFilename << std::endl;
  std::cout << "Propagation image: " << propagationFilename << std::endl;
  std::cout << "Patch radius: " << patchRadius << std::endl;
  std::cout << "Output: " << outputFilename << std::endl;

  typedef ITKHelpersTypes::UnsignedCharVectorImageType ImageType;
  ImageType::Pointer image = ImageType::New();
  ITKHelpers::ReadImage(imageFilename, image.GetPointer());

  Mask::Pointer mask = Mask::New();
  mask->Read(maskFilename);

  typedef StructurePropagation<ImageType>::PropagationLineImageType PropagationLineImageType;
  PropagationLineImageType::Pointer propagationLineImage = PropagationLineImageType::New();
  ITKHelpers::ReadImage(propagationFilename, propagationLineImage.GetPointer());

  StructurePropagation<ImageType> structurePropagation;
  structurePropagation.SetPropagationLineImage(propagationLineImage);
  structurePropagation.SetImage(image);
  structurePropagation.SetMask(mask);
  structurePropagation.SetPatchRadius(patchRadius);
  //structurePropagation.ComputeRegionsPropagationThroughHole();
  structurePropagation.ComputeRegionsPropagationAroundHole();
  structurePropagation.PropagateStructure();

  ITKHelpers::WriteImage(structurePropagation.GetOutputImage(), outputFilename);

  return EXIT_SUCCESS;
}

