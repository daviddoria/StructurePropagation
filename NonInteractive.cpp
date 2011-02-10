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

#include "Types.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"

int main(int argc, char* argv[])
{
  if(argc != 4)
    {
    std::cerr << "Required: image mask tileRadius" << std::endl;
    exit(-1);
    }

  typedef itk::ImageFileReader<ColorImageType> ImageReaderType;
  typedef itk::ImageFileReader<MaskImageType> MaskReaderType;

  std::string imageFilename = argv[1];
  std::cout << "Reading " << imageFilename << std::endl;

  ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader->SetFileName(imageFilename);
  imageReader->Update();

  std::string maskFilename = argv[2];
  std::cout << "Reading " << maskFilename << std::endl;

  std::string strRadius = argv[3];
  std::stringstream ss;
  ss << strRadius;
  int radius;
  ss >> radius;

  MaskReaderType::Pointer maskReader = MaskReaderType::New();
  maskReader->SetFileName(maskFilename);
  maskReader->Update();

  return EXIT_SUCCESS;
}

