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

#include "Helpers.h"
#include "Types.h"

// STL
#include <set>

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

std::vector<itk::Index<2> > FindIntersections(UnsignedCharScalarImageType::Pointer image)
{
  // We defined an intersection as a pixel which has more than 2 non-zero neighbors
  std::vector<itk::Index<2> > intersections;

  itk::Size<2> radius;
  radius.Fill(1);

  typedef itk::ConstNeighborhoodIterator<UnsignedCharScalarImageType> NeighborhoodIteratorType;
  NeighborhoodIteratorType iterator(radius, image,
                                    image->GetLargestPossibleRegion());

  std::vector<NeighborhoodIteratorType::OffsetType> neighbors = Helpers::CreateNeighborOffsetsInRadius(radius);

  while(!iterator.IsAtEnd())
    {
    if(iterator.GetCenterPixel() <= 0)
      {
      ++iterator;
      continue;
      }
    unsigned int neighborCounter = 0;
    for(unsigned int i = 0; i < neighbors.size(); i++)
      {
      bool inBounds;
      iterator.GetPixel(neighbors[i], inBounds);
      if(inBounds)
        {
        if(iterator.GetPixel(neighbors[i]) > 0)
          {
          neighborCounter++;
          }
        }
      }
    if(neighborCounter > 2)
      {
      std::cout << neighborCounter << " neighbors." << std::endl;
      intersections.push_back(iterator.GetIndex());
      }
    ++iterator;
    }
  return intersections;
}

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


unsigned int CountNumberOfNonNegativeValues(IntScalarImageType::Pointer image)
{
  itk::ImageRegionConstIterator<IntScalarImageType> imageIterator(image, image->GetLargestPossibleRegion());

  std::set<int> values;
  while(!imageIterator.IsAtEnd())
    {
    if(imageIterator.Get() >= 0)
      {
      values.insert(imageIterator.Get());
      }
    ++imageIterator;
    }

  return values.size();
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

void IndicesToBinaryImage(std::vector<itk::Index<2> > indices, itk::ImageRegion<2> imageRegion, UnsignedCharScalarImageType::Pointer image)
{
  //std::cout << "Setting " << indices.size() << " points to non-zero." << std::endl;

  // Blank the image
  image->SetRegions(imageRegion);
  image->Allocate();
  image->FillBuffer(0);

  // Set the pixels of indices in list to 255
  for(unsigned int i = 0; i < indices.size(); i++)
    {
    image->SetPixel(indices[i], 255);
    }
}

void IndicesToBinaryImage(std::set<itk::Index<2>, IndexComparison > indices, itk::ImageRegion<2> imageRegion,
                          UnsignedCharScalarImageType::Pointer image)
{
  // Blank the image
  image->SetRegions(imageRegion);
  image->Allocate();
  image->FillBuffer(0);

  // Set the pixels of indices in list to 255
  for(std::set<itk::Index<2>, IndexComparison >::iterator iterator = indices.begin(); iterator != indices.end(); iterator++)
    {
    image->SetPixel(*iterator, 255);
    }
}


void MaskImage(vtkSmartPointer<vtkImageData> VTKImage, vtkSmartPointer<vtkImageData> VTKSegmentMask, vtkSmartPointer<vtkImageData> VTKMaskedImage)
{
  int* dims = VTKImage->GetDimensions();

  VTKMaskedImage->SetDimensions(dims);
  VTKMaskedImage->AllocateScalars(VTK_UNSIGNED_CHAR, 4); // we might not have wanted to allocate here

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

  outputImage->SetDimensions(image->GetLargestPossibleRegion().GetSize()[0],
                             image->GetLargestPossibleRegion().GetSize()[1],
                             1);

  outputImage->AllocateScalars(VTK_UNSIGNED_CHAR, 1);

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

void ApplyTransparencyMask(vtkImageData* image, UnsignedCharScalarImageType::Pointer mask, vtkImageData* maskedImage)
{
  int dims[3];
  image->GetDimensions(dims);

  maskedImage->SetDimensions(dims);
  maskedImage->AllocateScalars(VTK_UNSIGNED_CHAR, 4);

  for(unsigned int i = 0; i < static_cast<unsigned int>(dims[0]); i++)
    {
    for(unsigned int j = 0; j < static_cast<unsigned int>(dims[1]); j++)
      {
      itk::Index<2> pixel;
      pixel[0] = i;
      pixel[1] = j;

      unsigned char* value = static_cast<unsigned char*>(image->GetScalarPointer(pixel[0], pixel[1],0));

      if(mask->GetPixel(pixel)) // if the pixel is in the masked region
        {
        value[3] = 255; // we want it to be visible (opaque)
        }
      else
        {
        value[3] = 0; // we want it to be invisible (transparent)
        }
      } // end j loop
    } // end i loop
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

std::vector<itk::Offset<2> > CreateNeighborOffsetsInRadius(itk::Size<2> radius)
{
  std::vector<itk::Offset<2> > neighbors;

  for(int i = -1 * static_cast<int>(radius[0]); i <= static_cast<int>(radius[0]); i++)
    {
    for(int j = -1 * static_cast<int>(radius[1]); j <= static_cast<int>(radius[1]); j++)
      {
      if(!(i==0 && j==0))
        {
        itk::Offset<2> offset = {{i,j}};
        neighbors.push_back(offset);
        } // end if
      } // end j loop
    } // end i loop
  return neighbors;
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


void GetDilatedImage(UnsignedCharScalarImageType::Pointer inputImage, UnsignedCharScalarImageType::Pointer dilatedImage,
                     unsigned int dilationRadius)
{
  typedef itk::BinaryBallStructuringElement<UnsignedCharScalarImageType::PixelType,2> StructuringElementType;
  StructuringElementType structuringElement;
  structuringElement.SetRadius(dilationRadius);
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
  numberedPixels->SetRegions(binaryImage->GetLargestPossibleRegion());
  numberedPixels->Allocate();
  numberedPixels->FillBuffer(-1);

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

void GrowNumberedPixelClusters(IntScalarImageType::Pointer numberedPixels, IntScalarImageType::Pointer clusteredPixels, unsigned int numberToSkip)
{
  std::cout << "Starting GrowNumberedPixelClusters..." << std::endl;

  // Input: numberedPixels - an image where each valid pixel is an id >= 0
  // Output: clusteredPixels - a modified version of numberedPixels where the regions of each id has grown
  // Details: The function expand the region associated with each label


  // Initialize by making the output the same as the input
  DeepCopy<IntScalarImageType>(numberedPixels, clusteredPixels);

  std::set<int> valuesSet = GetNonNegativeUniqueValues<IntScalarImageType>(clusteredPixels);
  std::vector<unsigned int> valuesVector(valuesSet.begin(), valuesSet.end());

  for(unsigned int i = 0; i < valuesVector.size(); i++)
    {
    if(i % numberToSkip == 0)
      {
      continue; // don't blank the label of every 'numberToSkip' pixels
      }
    ChangeValue<IntScalarImageType>(clusteredPixels,valuesVector[i], -2); // mark pixels to be filled with -2
    }

  std::vector<itk::Index<2> > pixelsToFill = FindPixelsWithValue<IntScalarImageType>(clusteredPixels, -2);
  std::vector<itk::Index<2> > validPixels = GetNonNegativePixels(clusteredPixels);

  for(unsigned int currentPixelId = 0; currentPixelId < pixelsToFill.size(); currentPixelId++) // loop over line pixels
    {
    unsigned int closestIndexId = FindClosestIndex(validPixels, pixelsToFill[currentPixelId]);
    clusteredPixels->SetPixel(pixelsToFill[currentPixelId], clusteredPixels->GetPixel(validPixels[closestIndexId]));
    }

  RelabelSequential(clusteredPixels);

  std::cout << "Finished ClusterNumberedPixels!" << std::endl;
}




void GrowNumberedPixelClustersWithIntersections(IntScalarImageType::Pointer numberedPixels, IntScalarImageType::Pointer clusteredPixels,
                                                unsigned int numberToSkip, std::set<itk::Index<2>, IndexComparison > intersections)
{
  std::cout << "Starting GrowNumberedPixelClustersWithIntersections..." << std::endl;

  // Input: numberedPixels - an image where each valid pixel is an id >= 0
  // Output: clusteredPixels - a modified version of numberedPixels where the regions of each id has grown
  // Details: The function expand the region associated with each label

  // Initialize by making the output the same as the input
  DeepCopy<IntScalarImageType>(numberedPixels, clusteredPixels);

  std::set<int> valuesSet = GetNonNegativeUniqueValues<IntScalarImageType>(clusteredPixels);
  std::vector<unsigned int> valuesVector(valuesSet.begin(), valuesSet.end());

  // Change labels ignoring intersections
  for(unsigned int i = 0; i < valuesVector.size(); i++)
    {
    if(i % numberToSkip == 0)
      {
      continue; // don't blank the label of every 'numberToSkip' pixels
      }
    ChangeValue<IntScalarImageType>(clusteredPixels,valuesVector[i], -2); // mark pixels to be filled with -2
    }

  // Change labels of intersections back to the original labels
  for(std::set<itk::Index<2>, IndexComparison>::iterator iterator = intersections.begin(); iterator != intersections.end(); iterator++)
    {
    clusteredPixels->SetPixel(*iterator, numberedPixels->GetPixel(*iterator));
    }

  std::vector<itk::Index<2> > pixelsToFill = FindPixelsWithValue<IntScalarImageType>(clusteredPixels, -2);
  std::vector<itk::Index<2> > validPixels = GetNonNegativePixels(clusteredPixels);

  std::cout << "Filling " << pixelsToFill.size() << " pixels." << std::endl;
  for(unsigned int currentPixelId = 0; currentPixelId < pixelsToFill.size(); currentPixelId++) // loop over line pixels
    {
    //unsigned int closestIndexId = FindClosestIndex(validPixels, pixelsToFill[currentPixelId]);
    unsigned int closestIndexId = FindClosestIndexWithMaxGeoDistance<IntScalarImageType>(validPixels, pixelsToFill[currentPixelId], numberedPixels, numberToSkip*2);

    clusteredPixels->SetPixel(pixelsToFill[currentPixelId], clusteredPixels->GetPixel(validPixels[closestIndexId]));
    }

  RelabelSequential(clusteredPixels);

  std::cout << "Finished GrowNumberedPixelClustersWithIntersections!" << std::endl;
}


std::vector<itk::Index<2> > GetLabelCenters(IntScalarImageType::Pointer image)
{
  int maxValue = GetMaxValue<IntScalarImageType>(image);
  //std::cout << "GetLabelCenters: MaxValue = " << maxValue << std::endl;

  itk::ImageRegionIterator<IntScalarImageType> imageIterator(image, image->GetLargestPossibleRegion());
  std::vector<itk::Index<2> > labelCenters;

  for(int i = 0; i <= maxValue; i++)
    {
    std::vector<itk::Index<2> > labelPixels = FindPixelsWithValue<IntScalarImageType>(image, i);

    //std::cout << "There are " << labelPixels.size() << " pixels with label " << i << std::endl;

    labelCenters.push_back(AverageIndex(labelPixels));
    }

  return labelCenters;
}


std::vector<itk::Index<2> > GetLabelCentersWithIntersections(IntScalarImageType::Pointer image, std::vector<itk::Index<2> > intersections)
{
  int maxValue = GetMaxValue<IntScalarImageType>(image);
  //std::cout << "GetLabelCenters: MaxValue = " << maxValue << std::endl;
  std::cout << "GetLabelCenters: number of intersections = " << intersections.size() << std::endl;

  itk::ImageRegionIterator<IntScalarImageType> imageIterator(image, image->GetLargestPossibleRegion());
  std::vector<itk::Index<2> > labelCenters;

  for(int i = 0; i <= maxValue; i++)
    {
    std::vector<itk::Index<2> > labelPixels = FindPixelsWithValue<IntScalarImageType>(image, i);

    //std::cout << "There are " << labelPixels.size() << " pixels with label " << i << std::endl;
    bool containsIntersection = false;
    itk::Index<2> containedIntersection;
    for(unsigned int intersectionId = 0; intersectionId < intersections.size(); intersectionId++)
      {
      if(ContainsElement<itk::Index<2> >(labelPixels, intersections[intersectionId]))
        {
        std::cout << "Contains intersection!!!" << std::endl;
        containsIntersection = true;
        containedIntersection = intersections[intersectionId];
        break;
        }
      }
    if(containsIntersection)
      {
      labelCenters.push_back(containedIntersection);
      }
    else
      {
      labelCenters.push_back(AverageIndex(labelPixels));
      }
    }

  return labelCenters;
}


std::vector<itk::Index<2> > GetLabelCentersWithIntersections(IntScalarImageType::Pointer image,
                                                             std::set<itk::Index<2>, IndexComparison > intersections)
{
  int maxValue = GetMaxValue<IntScalarImageType>(image);
  //std::cout << "GetLabelCenters: MaxValue = " << maxValue << std::endl;
  std::cout << "GetLabelCenters: number of intersections = " << intersections.size() << std::endl;

  itk::ImageRegionIterator<IntScalarImageType> imageIterator(image, image->GetLargestPossibleRegion());
  std::vector<itk::Index<2> > labelCenters;

  for(int i = 0; i <= maxValue; i++)
    {
    std::vector<itk::Index<2> > labelPixels = FindPixelsWithValue<IntScalarImageType>(image, i);

    //std::cout << "There are " << labelPixels.size() << " pixels with label " << i << std::endl;
    bool containsIntersection = false;
    itk::Index<2> containedIntersection;
    for(std::set<itk::Index<2>, IndexComparison>::iterator iterator = intersections.begin(); iterator != intersections.end(); iterator++)
      {
      if(ContainsElement<itk::Index<2> >(labelPixels, *iterator))
        {
        std::cout << "Contains intersection!!!" << std::endl;
        containsIntersection = true;
        containedIntersection = *iterator;
        break;
        }
      }
    if(containsIntersection)
      {
      labelCenters.push_back(containedIntersection);
      }
    else
      {
      labelCenters.push_back(AverageIndex(labelPixels));
      }
    }

  return labelCenters;
}

itk::Index<2> AverageIndex(std::vector<itk::Index<2> > pixels)
{
  //std::cout << "There are " << pixels.size() << " to average." << std::endl;
  int x = 0;
  int y = 0;
  for(unsigned int i = 0; i < pixels.size(); i++)
    {
    x += pixels[i][0];
    y += pixels[i][1];
    }
  itk::Index<2> center;
  center[0] = x/pixels.size();
  center[1] = y/pixels.size();

  return center;
}

void RelabelSequential(IntScalarImageType::Pointer image)
{
  //std::cout << "Starting RelabelSequential..." << std::endl;
  // Add all of the values in the image to a set
  itk::ImageRegionConstIterator<IntScalarImageType> inputIterator(image, image->GetLargestPossibleRegion());

  std::set<int> values;
  while(!inputIterator.IsAtEnd())
    {
    if(inputIterator.Get() >= 0)
      {
      values.insert(inputIterator.Get());
      }
    ++inputIterator;
    }

  itk::Image<bool, 2>::Pointer visitedImage = itk::Image<bool, 2>::New();
  visitedImage->SetRegions(image->GetLargestPossibleRegion());
  visitedImage->Allocate();
  visitedImage->FillBuffer(0);

  // Relabel each of the pixels with matching labels starting at 0
  int labelId = 0;
  for(std::set<int>::iterator it = values.begin(); it != values.end(); it++)
    {
    //std::cout << "Looking for pixels with value " << *it << std::endl;
    itk::ImageRegionIterator<IntScalarImageType> outputIterator(image, image->GetLargestPossibleRegion());
    outputIterator.GoToBegin();
    unsigned int counter = 0;
    while(!outputIterator.IsAtEnd())
      {
      if(outputIterator.Get() == *it)
        {
        if(!visitedImage->GetPixel(outputIterator.GetIndex()))
          {
          //std::cout << outputIterator.Get() << " equals " << *it << std::endl;
          outputIterator.Set(labelId);
          visitedImage->SetPixel(outputIterator.GetIndex(), 1); // mark as visited
          counter++;
          }
        }
      ++outputIterator;
      }
    //std::cout << "RelabelSequential: Marked " << counter << " pixels with label " << labelId << std::endl;
    if(counter == 0)
      {
      std::cerr << "Every label should have at least one pixel!" << std::endl;
      exit(-1);
      }
    //std::cout << "Changed " << counter << " pixels from " << *it << " to " << labelId << std::endl;
    labelId++;
    }

  // For testing only
/*
  for(unsigned int i = 0; i < values.size(); i++)
    {
    unsigned int number = CountPixelsWithValue<IntScalarImageType>(image, i);
    std::cout << "RelabelSequential Test: Label " << i << " has " << number << " pixels." << std::endl;
    }
*/
  //std::cout << "Finished RelabelSequential!" << std::endl;
}

void WritePixels(itk::Size<2> imageSize, std::vector<itk::Index<2> > pixels, std::string filename)
{
  UnsignedCharScalarImageType::Pointer image = UnsignedCharScalarImageType::New();
  itk::Index<2> start;
  start.Fill(0);
  itk::ImageRegion<2> region(start, imageSize);
  image->SetRegions(region);
  image->Allocate();
  image->FillBuffer(0);

  //std::cout << "WritePixels: Writing " << pixels.size() << " pixels." << std::endl;

  for(unsigned int i = 0; i < pixels.size(); i++)
    {
    image->SetPixel(pixels[i], 255);
    }

  WriteImage<UnsignedCharScalarImageType>(image, filename);
}

double DistanceBetweenIndices(itk::Index<2> pixel1, itk::Index<2> pixel2)
{
  itk::Point<double,2> p1;
  p1[0] = pixel1[0];
  p1[1] = pixel1[1];

  itk::Point<double,2> p2;
  p2[0] = pixel2[0];
  p2[1] = pixel2[1];

  return p1.EuclideanDistanceTo(p2);
}

unsigned int FindClosestIndex(std::vector<itk::Index<2> > listOfPixels, itk::Index<2> queryPixel)
{
  double minimumDistance = itk::NumericTraits<double>::max();
  unsigned int closestIndex = 0;
  for(unsigned int i = 0; i < listOfPixels.size(); i++)
    {
    double distance = DistanceBetweenIndices(listOfPixels[i], queryPixel);
    if(distance < minimumDistance)
      {
      minimumDistance = distance;
      closestIndex = i;
      }
    }
  return closestIndex;
}


int GetFirstNonNegativeValidNeighbor(IntScalarImageType::Pointer image, itk::Index<2> pixel)
{

  itk::Size<2> radius;
  radius.Fill(1);

  // Create a single pixel region - we only want to traverse the neighborhood of the iterator once
  itk::ImageRegion<2> region(pixel, radius);

  typedef itk::NeighborhoodIterator<IntScalarImageType> NeighborhoodIteratorType;
  NeighborhoodIteratorType iterator(radius, image, region);

  std::vector<NeighborhoodIteratorType::OffsetType> neighbors = Helpers::CreateNeighborOffsetsInRadius(radius);

  while(!iterator.IsAtEnd()) // this should only be one time
    {
    for(unsigned int i = 0; i < neighbors.size(); i++)
      {
      bool inBounds;
      iterator.GetPixel(neighbors[i], inBounds);
      if(inBounds)
        {
        std::cout << "Pixel " << i << " " << iterator.GetPixel(neighbors[i]) << std::endl;
        if(iterator.GetPixel(neighbors[i]) >= 0)
          {
          return iterator.GetPixel(neighbors[i]);
          }
        }
      }
    ++iterator;
    }

  std::cerr << "GetFirstNonNegativeValidNeighbor: No valid neighbors!" << std::endl;
  exit(-1);
}

std::vector<itk::Index<2> > GetNonNegativeValidNeighbors(IntScalarImageType::Pointer image, itk::Index<2> pixel)
{
  //std::cout << "Starting GetNonNegativeValidNeighbors..." << std::endl;
  std::vector<itk::Index<2> > nonNegativeValidNeighbors;

  if(image->GetPixel(pixel) < 0)
    {
    std::cout << "GetNonNegativeValidNeighbors: Zero neighbors!" << std::endl;
    return nonNegativeValidNeighbors;
    }

  itk::Size<2> radius;
  radius.Fill(1);

  // Create a single pixel region - we only want to traverse the neighborhood of the iterator once
  itk::ImageRegion<2> region(pixel, radius);

  typedef itk::NeighborhoodIterator<IntScalarImageType> NeighborhoodIteratorType;
  NeighborhoodIteratorType iterator(radius, image, region);

  std::vector<NeighborhoodIteratorType::OffsetType> neighbors = Helpers::CreateNeighborOffsetsInRadius(radius);

  while(!iterator.IsAtEnd()) // this should only be one time
    {
    for(unsigned int i = 0; i < neighbors.size(); i++)
      {
      bool inBounds;
      iterator.GetPixel(neighbors[i], inBounds);
      if(inBounds)
        {
        if(iterator.GetPixel(neighbors[i]) >= 0)
          {
          nonNegativeValidNeighbors.push_back(iterator.GetIndex() + neighbors[i]);
          }
        }
      }
    ++iterator;
    }

  //std::cout << "Finished GetNonNegativeValidNeighbors!" << std::endl;

  return nonNegativeValidNeighbors;
}

std::vector<itk::Index<2> > GetNonNegativePixels(IntScalarImageType::Pointer image)
{
  std::vector<itk::Index<2> > nonNegativePixels;

  itk::ImageRegionConstIterator<IntScalarImageType> imageIterator(image, image->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    if(imageIterator.Get() >= 0)
      {
      nonNegativePixels.push_back(imageIterator.GetIndex());
      }

    ++imageIterator;
    }
  return nonNegativePixels;
}

} // end namespace

bool operator<(IndexDistance object1, IndexDistance object2)
{
  return object1.Distance < object2.Distance;
}
