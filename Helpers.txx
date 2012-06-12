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
#include "itkCastImageFilter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageFileWriter.h"
#include "itkNeighborhoodIterator.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkScalarToRGBColormapImageFilter.h"

// Boost
#include <boost/graph/graph_concepts.hpp>

// Custom
#include "DijkstraBinaryImage.h"

namespace Helpers
{

template <typename TImage>
typename TImage::PixelType GetMaxValue(typename TImage::Pointer image)
{
  typedef itk::MinimumMaximumImageCalculator<TImage>
          ImageCalculatorFilterType;

  typename ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New ();
  imageCalculatorFilter->SetImage(image);
  imageCalculatorFilter->Compute();
  return imageCalculatorFilter->GetMaximum();
}

template <typename TImage>
unsigned int CountPixelsWithValue(typename TImage::Pointer image, typename TImage::PixelType pixel)
{
  itk::ImageRegionConstIterator<TImage> imageIterator(image, image->GetLargestPossibleRegion());
  unsigned int counter = 0;
  while(!imageIterator.IsAtEnd())
    {
    if(imageIterator.Get() == pixel)
      {
      counter++;
      }
    ++imageIterator;
    }
  return counter;
}


template <typename TImage>
std::vector<itk::Index<2> > FindPixelsWithValue(typename TImage::Pointer image, typename TImage::PixelType pixel)
{
  std::vector<itk::Index<2> > pixels;
  itk::ImageRegionConstIterator<TImage> imageIterator(image, image->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    if(imageIterator.Get() == pixel)
      {
      pixels.push_back(imageIterator.GetIndex());
      }
    ++imageIterator;
    }
  return pixels;
}

template <typename TImage>
void DeepCopy(typename TImage::Pointer input, typename TImage::Pointer output)
{
  itk::ImageRegionConstIterator<TImage> inputIterator(input, input->GetLargestPossibleRegion());
  output->SetRegions(input->GetLargestPossibleRegion());
  output->Allocate();
  itk::ImageRegionIterator<TImage> outputIterator(output, output->GetLargestPossibleRegion());

  while(!inputIterator.IsAtEnd())
    {
    outputIterator.Set(inputIterator.Get());
    ++inputIterator;
    ++outputIterator;
    }
}

// Convert an ITK image to a VTK image for display
template <typename TImageType>
void ITKImageToVTKImage(typename TImageType::Pointer image, vtkSmartPointer<vtkImageData> outputImage)
{
  // Setup and allocate the image data
  outputImage->SetNumberOfScalarComponents(TImageType::PixelType::GetNumberOfComponents());
  outputImage->SetScalarTypeToUnsignedChar();
  outputImage->SetDimensions(image->GetLargestPossibleRegion().GetSize()[0],
                             image->GetLargestPossibleRegion().GetSize()[1],
                             1);

  outputImage->AllocateScalars();

  // Copy all of the input image pixels to the output image
  itk::ImageRegionConstIteratorWithIndex<TImageType> imageIterator(image,image->GetLargestPossibleRegion());
  imageIterator.GoToBegin();

  while(!imageIterator.IsAtEnd())
    {
    unsigned char* pixel = static_cast<unsigned char*>(outputImage->GetScalarPointer(imageIterator.GetIndex()[0],
                                                                                     imageIterator.GetIndex()[1],0));
    for(unsigned int component = 0; component < TImageType::PixelType::GetNumberOfComponents(); component++)
      {
      pixel[component] = static_cast<unsigned char>(imageIterator.Get()[component]);
      }

    ++imageIterator;
    }
}

template <typename TImageType>
void WriteImage(typename TImageType::Pointer image, std::string filename)
{
  typedef  itk::ImageFileWriter< TImageType > WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(filename);
  writer->SetInput(image);
  writer->Update();
}

template <typename TImage>
void WriteColorMappedImage(typename TImage::Pointer image, std::string filename)
{
  typedef itk::RescaleIntensityImageFilter< TImage, UnsignedCharScalarImageType > RescaleFilterType;
  typename RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  rescaleFilter->SetInput(image);
  rescaleFilter->SetOutputMinimum(0);
  rescaleFilter->SetOutputMaximum(255);
  rescaleFilter->Update();

  typedef itk::RGBPixel<unsigned char>    RGBPixelType;
  typedef itk::Image<RGBPixelType, 2>  RGBImageType;

  typedef itk::ScalarToRGBColormapImageFilter<UnsignedCharScalarImageType, RGBImageType> RGBFilterType;
  RGBFilterType::Pointer rgbfilter = RGBFilterType::New();
  rgbfilter->SetInput(rescaleFilter->GetOutput());
  //rgbfilter->SetColormap( RGBFilterType::Hot );
  rgbfilter->SetColormap( RGBFilterType::Jet );
  rgbfilter->Update();

  typedef  itk::ImageFileWriter<RGBImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(filename);
  writer->SetInput(rgbfilter->GetOutput());
  writer->Update();
}


template<typename TImage>
void CastAndWriteImage(typename TImage::Pointer image, std::string filename)
{
  typedef itk::CovariantVector<unsigned char,TImage::PixelType::Dimension> UnsignedCharVectorPixelType;
  typedef itk::Image<UnsignedCharVectorPixelType, 2> UnsignedCharVectorImageType;

  typedef itk::CastImageFilter< TImage, UnsignedCharVectorImageType > CastFilterType;
  typename CastFilterType::Pointer castFilter = CastFilterType::New();
  castFilter->SetInput(image);
  castFilter->Update();

  WriteImage<UnsignedCharVectorImageType>(castFilter->GetOutput(), filename);
}

template <typename TPixelType>
double SquaredDifference(TPixelType pixel1, TPixelType pixel2)
{
  return (pixel1-pixel2).GetSquaredNorm();
}

template <typename TImage>
void CopyPatchIntoImage(typename TImage::Pointer sourceImage, typename TImage::Pointer targetImage,
                                   itk::ImageRegion<2> sourceRegion, itk::ImageRegion<2> targetRegion)
{
  itk::ImageRegionConstIterator<TImage> sourceImageIterator(sourceImage, sourceRegion);
  itk::ImageRegionIterator<TImage> targetImageIterator(targetImage, targetRegion);

  while(!sourceImageIterator.IsAtEnd())
    {
    targetImageIterator.Set(sourceImageIterator.Get());
    ++sourceImageIterator;
    ++targetImageIterator;
    }
}

template <typename TImage>
void CopySelfPatchIntoTargetRegion(typename TImage::Pointer image, MaskImageType::Pointer mask,
                                   itk::ImageRegion<2> sourceRegion, itk::ImageRegion<2> targetRegion)
{
  //assert(IsCompletelyValid(sourceRegion));

  itk::ImageRegionConstIterator<TImage> imageSourceIterator(image, sourceRegion);
  itk::ImageRegionIterator<TImage> imageTargetIterator(image, targetRegion);
  itk::ImageRegionConstIterator<MaskImageType> maskIterator(mask, targetRegion);

  while(!imageSourceIterator.IsAtEnd())
    {
    if(maskIterator.Get()) // we are in the target region
      {
      imageTargetIterator.Set(imageSourceIterator.Get());
      }
    ++imageSourceIterator;
    ++imageTargetIterator;
    ++maskIterator;
    }
}

template <typename TImage>
void ChangeValue(typename TImage::Pointer image, typename TImage::PixelType sourceValue, typename TImage::PixelType targetValue)
{
  itk::ImageRegionIterator<TImage> imageIterator(image, image->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    if(imageIterator.Get() == sourceValue)
      {
      imageIterator.Set(targetValue);
      }
    ++imageIterator;
    }
}

template <typename TImage>
std::set<typename TImage::PixelType> GetUniqueValues(typename TImage::Pointer image)
{
  std::set<typename TImage::PixelType> values;
  itk::ImageRegionConstIterator<TImage> imageIterator(image, image->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    values.insert(imageIterator.Get());
    ++imageIterator;
    }
  return values;
}


template <typename TImage>
std::set<typename TImage::PixelType> GetNonNegativeUniqueValues(typename TImage::Pointer image)
{
  std::set<typename TImage::PixelType> values;
  itk::ImageRegionConstIterator<TImage> imageIterator(image, image->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    if(imageIterator.Get() >= 0)
      {
      values.insert(imageIterator.Get());
      }
    ++imageIterator;
    }
  return values;
}

template <typename T>
bool ContainsElement(std::vector<T> vec, T element)
{
  for(unsigned int i = 0; i < vec.size(); i++)
    {
    if(vec[i] == element)
      {
      return true;
      }
    }
  return false;
}

template <typename T>
bool HasNeighborWithValue(typename T::Pointer image, itk::Index<2> pixel, typename T::PixelType value)
{
  itk::Size<2> radius;
  radius.Fill(1);

  // Create a single pixel region - we only want to traverse the neighborhood of the iterator once
  itk::ImageRegion<2> region(pixel, radius);

  typedef itk::NeighborhoodIterator<T> NeighborhoodIteratorType;
  NeighborhoodIteratorType iterator(radius, image, region);

  //std::vector<NeighborhoodIteratorType::OffsetType> neighbors = CreateNeighborOffsetsInRadius(radius);
  std::vector<itk::Offset<2> > neighbors = CreateNeighborOffsetsInRadius(radius);

  while(!iterator.IsAtEnd()) // this should only be one time
    {
    for(unsigned int i = 0; i < neighbors.size(); i++)
      {
      bool inBounds;
      iterator.GetPixel(neighbors[i], inBounds);
      if(inBounds)
        {
        if(iterator.GetPixel(neighbors[i]) == value)
          {
          return true;
          }
        }
      }
    ++iterator;
    }
  return false;
}

template <typename T>
unsigned int FindClosestIndexGeodesic(typename T::Pointer image, std::vector<itk::Index<2> > listOfPixels,
                                      itk::Index<2> queryPixel)
{
  std::cout << "There are " << listOfPixels.size() << " pixels to compute the path to." << std::endl;

  int minimumDistance = itk::NumericTraits<int>::max();
  unsigned int closestIndex = 0;
  for(unsigned int i = 0; i < listOfPixels.size(); i++)
    {
    DijkstraBinaryImage<T> DistanceFilter;
    DistanceFilter.SetImage(image);

    itk::Index<2> startPoint = queryPixel;
    itk::Index<2> endPoint = listOfPixels[i];

    int distance = DistanceFilter.ComputeShortestPath(startPoint, endPoint).size();

    if(distance < minimumDistance)
      {
      minimumDistance = distance;
      closestIndex = i;
      }
    }
  return closestIndex;
}

template <typename T>
unsigned int FindClosestIndexWithMaxGeoDistance(std::vector<itk::Index<2> > listOfPixels, itk::Index<2> queryPixel,
                                                typename T::Pointer image, int maxDistance)
{
  std::vector<IndexDistance> allDistances;

  for(unsigned int i = 0; i < listOfPixels.size(); i++)
    {
    IndexDistance indexDistance;
    double distance = DistanceBetweenIndices(listOfPixels[i], queryPixel);
    indexDistance.Distance = distance;
    indexDistance.Index = i;
    allDistances.push_back(indexDistance);
    }

  std::sort(allDistances.begin(), allDistances.end());

  for(int i = 0; i < allDistances.size(); i++)
    {
    DijkstraBinaryImage<T> dijkstra;
    dijkstra.SetImage(image);
    itk::Index<2> startPoint = queryPixel;
    itk::Index<2> endPoint = listOfPixels[allDistances[i].Index];

    int distance = dijkstra.ComputeShortestPath(startPoint, endPoint).size();
    if(distance < maxDistance)
      {
      return allDistances[i].Index;
      }
    }

  // Should never get here!
  std::cerr << "FindClosestIndexWithMaxGeoDistance: There were no pixels with distance < maxDistance (" << maxDistance << ")!" << std::endl;
  std::cout << "FindClosestIndexWithMaxGeoDistance: Closest is " << allDistances[0].Index << " at distance " << allDistances[0].Distance << std::endl;
  //exit(-1);
  return 0;
}


} // end namespace