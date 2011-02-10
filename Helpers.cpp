#include "Helpers.h"
#include "Types.h"

// ITK
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBresenhamLine.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileWriter.h"

// VTK
#include <vtkPolyData.h>

namespace Helpers
{

bool IsIntersectTargetRegion(itk::ImageRegion<2> patch, UnsignedCharScalarImageType::Pointer mask)
{
  itk::ImageRegionConstIterator<UnsignedCharScalarImageType> imageIterator(mask, patch);

  while(!imageIterator.IsAtEnd())
    {
    if(imageIterator.Get() > 0) // the mask is non-zero
      {
      return true;
      }
    ++imageIterator;
    }

  return false;
}

bool IsValidPatch(itk::ImageRegion<2> patch, UnsignedCharScalarImageType::Pointer mask)
{
  // The patch must be entirely within the image
  if(!mask->GetLargestPossibleRegion().IsInside(patch))
    {
    return false;
    }

  // There must be zero masked (invalid) pixels
  itk::ImageRegionConstIterator<UnsignedCharScalarImageType> imageIterator(mask, patch);

  bool validRegion = true;
  while(!imageIterator.IsAtEnd())
    {
    if(imageIterator.Get() > 0) // the mask is non-zero
      {
      return false;
      }
    ++imageIterator;
    }

  return true;
}

itk::Index<2> GetCorner(itk::Index<2> center, unsigned int radius)
{
  itk::Index<2> corner;
  corner[0] = center[0] - radius;
  corner[1] = center[1] - radius;

  return corner;
}

itk::ImageRegion<2> GetCenteredRegionRadius(itk::Index<2> centerPixel, unsigned int radius)
{
  itk::Size<2> size;
  size.Fill(radius*2 + 1);

  itk::Index<2> corner;
  corner[0] = centerPixel[0] - radius;
  corner[1] = centerPixel[1] - radius;

  return itk::ImageRegion<2>(corner, size);
}

double SumOfShortestDistances(std::vector<itk::Index<2> > indices1, std::vector<itk::Index<2> > indices2)
{
  // not finished
  for(unsigned int i = 0; i < indices1.size(); i++)
    {
    for(unsigned int j = 0; j < indices1.size(); j++)
      {
      //double dist = p0.EuclideanDistanceTo(p1);
      }
    }
  return 0;
}

double SumOfShortestDistances(std::vector<itk::Point<int, 2> > points1, std::vector<itk::Point<int, 2> > points2)
{
  double totalDistance = 0;
  for(unsigned int i = 0; i < points1.size(); i++)
    {
    double minDist = 1e8;
    double dist = 0;
    for(unsigned int j = 0; j < points1.size(); j++)
      {
      dist = points1[i].EuclideanDistanceTo(points2[j]);
      if(dist < minDist)
        {
        minDist = dist;
        }
      }
    totalDistance += dist;
    }
  return totalDistance;
}

void WriteWhitePatches(itk::ImageRegion<2> imageRegion, std::vector<itk::ImageRegion<2> > regions, std::string filename)
{
  std::cout << "Writing " << regions.size() << " patches." << std::endl;

  UnsignedCharScalarImageType::Pointer image = UnsignedCharScalarImageType::New();
  image->SetRegions(imageRegion);
  image->Allocate();
  image->FillBuffer(0);

  for(unsigned int i = 0; i < regions.size(); i++)
    {
    itk::ImageRegionIterator<UnsignedCharScalarImageType> imageIterator(image,regions[i]);

    while(!imageIterator.IsAtEnd())
      {
      imageIterator.Set(255);

      ++imageIterator;
      }

    }

  typedef  itk::ImageFileWriter< UnsignedCharScalarImageType  > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(filename);
  writer->SetInput(image);
  writer->Update();
}

void IndicesToBinaryImage(std::vector<itk::Index<2> > indices, UnsignedCharScalarImageType::Pointer image)
{
  //std::cout << "Setting " << indices.size() << " points to non-zero." << std::endl;

  // Blank the image
  image->FillBuffer(0);

  // Set the pixels of indices in list to 255
  for(unsigned int i = 0; i < indices.size(); i++)
    {
    image->SetPixel(indices[i], 255);
    }
}


void MaskImage(vtkSmartPointer<vtkImageData> VTKImage, vtkSmartPointer<vtkImageData> VTKSegmentMask, vtkSmartPointer<vtkImageData> VTKMaskedImage)
{
  int* dims = VTKImage->GetDimensions();

  VTKMaskedImage->SetDimensions(dims);
  VTKMaskedImage->SetNumberOfScalarComponents(4);
  VTKMaskedImage->SetScalarTypeToUnsignedChar();

  // int dims[3]; // can't do this
  for (int y = 0; y < dims[1]; y++)
    {
    for (int x = 0; x < dims[0]; x++)
      {

      unsigned char* imagePixel = static_cast<unsigned char*>(VTKImage->GetScalarPointer(x,y,0));
      unsigned char* maskPixel = static_cast<unsigned char*>(VTKSegmentMask->GetScalarPointer(x,y,0));
      unsigned char* outputPixel = static_cast<unsigned char*>(VTKMaskedImage->GetScalarPointer(x,y,0));

      outputPixel[0] = imagePixel[0];

      if(VTKImage->GetNumberOfScalarComponents() == 3)
        {
        outputPixel[1] = imagePixel[1];
        outputPixel[2] = imagePixel[2];
        }
      else // Grayscale should have all components equal to the first component
        {
        outputPixel[1] = imagePixel[0];
        outputPixel[2] = imagePixel[0];
        }

      if(maskPixel[0] == 0)
        {
        outputPixel[3] = 0;
        }
      else
        {
        outputPixel[3] = 255;
        }

      }
    }
}


// Specialization for image types with pixel types without [] operator
template <>
void ITKImageToVTKImage<UnsignedCharScalarImageType>(UnsignedCharScalarImageType::Pointer image, vtkSmartPointer<vtkImageData> outputImage)
{
  // Setup and allocate the image data
  outputImage->SetNumberOfScalarComponents(1);
  outputImage->SetScalarTypeToUnsignedChar();
  outputImage->SetDimensions(image->GetLargestPossibleRegion().GetSize()[0],
                             image->GetLargestPossibleRegion().GetSize()[1],
                             1);

  outputImage->AllocateScalars();

  // Copy all of the input image pixels to the output image
  itk::ImageRegionConstIteratorWithIndex<MaskImageType> imageIterator(image,image->GetLargestPossibleRegion());
  imageIterator.GoToBegin();

  while(!imageIterator.IsAtEnd())
    {
    unsigned char* pixel = static_cast<unsigned char*>(outputImage->GetScalarPointer(imageIterator.GetIndex()[0],
                                                                                     imageIterator.GetIndex()[1],0));
    pixel[0] = static_cast<unsigned char>(imageIterator.Get());
    ++imageIterator;
    }
}

void VTKImageToITKImage(vtkImageData* image, UnsignedCharScalarImageType::Pointer outputImage)
{
  if(image->GetNumberOfScalarComponents() > 1)
    {
    std::cerr << "VTKImagetoITKImage only handles images with 1 component!" << std::endl;
    return;
    }
  int dims[3];
  image->GetDimensions(dims);

  itk::Size<2> size;
  size[0] = dims[0];
  size[1] = dims[1];

  itk::Index<2> index;
  index.Fill(0);

  itk::ImageRegion<2> region(index, size);
  outputImage->SetRegions(region);
  outputImage->Allocate();

  for(unsigned int i = 0; i < static_cast<unsigned int>(dims[0]); i++)
    {
    for(unsigned int j = 0; j < static_cast<unsigned int>(dims[1]); j++)
      {
      itk::Index<2> pixel;
      pixel[0] = i;
      pixel[1] = j;

      unsigned char* value = static_cast<unsigned char*>(image->GetScalarPointer(pixel[0], pixel[1],0));
      outputImage->SetPixel(pixel, value[0]);
      }
    }
}


std::vector<itk::Index<2> > PolyDataToPixelList(vtkPolyData* polydata)
{
  // Convert vtkPoints to indices
  std::vector<itk::Index<2> > linePoints;
  for(vtkIdType i = 0; i < polydata->GetNumberOfPoints(); i++)
    {
    itk::Index<2> index;
    double p[3];
    polydata->GetPoint(i,p);
    index[0] = round(p[0]);
    index[1] = round(p[1]);
    linePoints.push_back(index);
    }

  // Compute the indices between every pair of points
  std::vector<itk::Index<2> > allIndices;
  for(unsigned int linePointId = 1; linePointId < linePoints.size(); linePointId++)
    {
    itk::Index<2> index0 = linePoints[linePointId-1];
    itk::Index<2> index1 = linePoints[linePointId];
    // Currently need the distance between the points for Bresenham (pending patch in Gerrit)
    itk::Point<float,2> point0;
    itk::Point<float,2> point1;
    for(unsigned int i = 0; i < 2; i++)
      {
      point0[i] = index0[i];
      point1[i] = index1[i];
      }
    float distance = point0.EuclideanDistanceTo(point1);
    itk::BresenhamLine<2> line;
    std::vector<itk::Offset<2> > offsets = line.BuildLine(point1-point0, distance);
    for(unsigned int i = 0; i < offsets.size(); i++)
      {
      allIndices.push_back(index0 + offsets[i]);
      }
    }

  return allIndices;
}


void DrawLineInImage(UnsignedCharScalarImageType::Pointer image, itk::Index<2> p0, itk::Index<2> p1)
{
  // Currently need the distance between the points for Bresenham (pending patch in Gerrit)
  itk::Point<float,2> point0;
  itk::Point<float,2> point1;
  for(unsigned int i = 0; i < 2; i++)
    {
    point0[i] = p0[i];
    point1[i] = p1[i];
    }
    
  float distance = point0.EuclideanDistanceTo(point1);
  itk::BresenhamLine<2> line;
  std::vector<itk::Offset<2> > offsets = line.BuildLine(point1-point0, distance);
  //std::cout << "There are " << offsets.size() << " pixels to turn white in the line." << std::endl;
  
  for(unsigned int i = 0; i < offsets.size(); i++)
    {
    image->SetPixel(p0 + offsets[i], 255);
    }

}

std::vector<itk::Index<2> > BinaryImageToPixelList(UnsignedCharScalarImageType::Pointer image)
{
  std::vector<itk::Index<2> > pixelList;

  itk::ImageRegionConstIterator<UnsignedCharScalarImageType> imageIterator(image, image->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    if(imageIterator.Get() > 0)
      {
      pixelList.push_back(imageIterator.GetIndex());
      }

    ++imageIterator;
    }
  return pixelList;
}

unsigned int CountNonZeroPixels(UnsignedCharScalarImageType::Pointer image)
{
  itk::ImageRegionConstIterator<UnsignedCharScalarImageType> imageIterator(image, image->GetLargestPossibleRegion());
  unsigned int count = 0;
  while(!imageIterator.IsAtEnd())
    {
    if(imageIterator.Get() > 0)
      {
      count++;
      }

    ++imageIterator;
    }

  return count;
}


void GetDilatedImage(UnsignedCharScalarImageType::Pointer inputImage, UnsignedCharScalarImageType::Pointer dilatedImage)
{
  typedef itk::BinaryBallStructuringElement<UnsignedCharScalarImageType::PixelType,2> StructuringElementType;
  StructuringElementType structuringElement;
  structuringElement.SetRadius(2);
  structuringElement.CreateStructuringElement();

  typedef itk::BinaryDilateImageFilter <UnsignedCharScalarImageType, UnsignedCharScalarImageType, StructuringElementType>
          BinaryDilateImageFilterType;

  BinaryDilateImageFilterType::Pointer dilateFilter = BinaryDilateImageFilterType::New();
  dilateFilter->SetInput(inputImage);
  dilateFilter->SetKernel(structuringElement);
  dilateFilter->Update();

  dilatedImage->Graft(dilateFilter->GetOutput());

}

void NumberPixels(UnsignedCharScalarImageType::Pointer binaryImage, IntScalarImageType::Pointer numberedPixels)
{
  itk::ImageRegionConstIterator<UnsignedCharScalarImageType> imageIterator(binaryImage, binaryImage->GetLargestPossibleRegion());
  itk::ImageRegionIterator<IntScalarImageType> numberedPixelsIterator(numberedPixels, numberedPixels->GetLargestPossibleRegion());

  int counter = 0;
  while(!imageIterator.IsAtEnd())
    {
    if(imageIterator.Get())
      {
      numberedPixelsIterator.Set(counter);
      counter++;
      }
    else
      {
      numberedPixelsIterator.Set(-1);
      }

    ++imageIterator;
    ++numberedPixelsIterator;
    }
}

} // end namespace
