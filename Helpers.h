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

#ifndef HELPERS_H
#define HELPERS_H

// Custom
#include "Types.h"
#include "IndexComparison.h"

// ITK
#include "itkImage.h"
#include "itkIndex.h"

// VTK
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkPolyData.h>
#include <vtkUnsignedCharArray.h>

// STL
#include <string>
#include <set>

namespace Helpers
{

////////////////// Non templated functions (defined in Helpers.cpp) //////////////////

void ApplyTransparencyMask(vtkImageData* image, UnsignedCharScalarImageType::Pointer mask, vtkImageData* maskedImage);

double DistanceBetweenIndices(itk::Index<2> pixel1, itk::Index<2> pixel2);

unsigned int FindClosestIndex(std::vector<itk::Index<2> > listOfPixels, itk::Index<2> queryPixel);

std::vector<itk::Index<2> > GetNonNegativeValidNeighbors(IntScalarImageType::Pointer image, itk::Index<2> pixel);

std::vector<itk::Index<2> > FindIntersections(UnsignedCharScalarImageType::Pointer image);

int GetFirstNonNegativeValidNeighbor(IntScalarImageType::Pointer image, itk::Index<2> pixel);

std::vector<itk::Index<2> > GetNonNegativePixels(IntScalarImageType::Pointer image);

unsigned int CountNumberOfNonNegativeValues(IntScalarImageType::Pointer image);

std::vector<itk::Index<2> > GetLabelCenters(IntScalarImageType::Pointer image);

std::vector<itk::Index<2> > GetLabelCentersWithIntersections(IntScalarImageType::Pointer image, std::vector<itk::Index<2> > intersections);

std::vector<itk::Index<2> > GetLabelCentersWithIntersections(IntScalarImageType::Pointer image,
                                                             std::set<itk::Index<2>, IndexComparison > intersections);

itk::Index<2> AverageIndex(std::vector<itk::Index<2> > pixels);

void DrawLineInImage(UnsignedCharScalarImageType::Pointer image, itk::Index<2> p0, itk::Index<2> p1);

void GrowNumberedPixelClusters(IntScalarImageType::Pointer numberedPixels, IntScalarImageType::Pointer clusteredPixels,
                               unsigned int numberToSkip);

void GrowNumberedPixelClustersWithIntersections(IntScalarImageType::Pointer numberedPixels, IntScalarImageType::Pointer clusteredPixels,
                               unsigned int numberToSkip, std::vector<itk::Index<2> > intersections);

void GrowNumberedPixelClustersWithIntersections(IntScalarImageType::Pointer numberedPixels, IntScalarImageType::Pointer clusteredPixels,
                               unsigned int numberToSkip, std::set<itk::Index<2>, IndexComparison> intersections);

void NumberPixels(UnsignedCharScalarImageType::Pointer binaryImage, IntScalarImageType::Pointer numberedPixels);

void RelabelSequential(IntScalarImageType::Pointer numberedPixels);

itk::Index<2> GetCorner(itk::Index<2> center, unsigned int radius);

std::vector<itk::Offset<2> > CreateNeighborOffsetsInRadius(itk::Size<2> radius);

void GetDilatedImage(UnsignedCharScalarImageType::Pointer inputImage, UnsignedCharScalarImageType::Pointer dilatedImage,
                     unsigned int dilationRadius);

double SumOfShortestDistances(std::vector<itk::Index<2> > indices1, std::vector<itk::Index<2> > indices2);

double SumOfShortestDistances(std::vector<itk::Point<int, 2> > points1, std::vector<itk::Point<int, 2> > points2);

// Mark each pixel at the specified 'indices' as a non-zero pixel in 'image'
void IndicesToBinaryImage(std::vector<itk::Index<2> > indices, itk::ImageRegion<2> imageRegion, UnsignedCharScalarImageType::Pointer image);
void IndicesToBinaryImage(std::set<itk::Index<2>, IndexComparison > indices, itk::ImageRegion<2> imageRegion,
                          UnsignedCharScalarImageType::Pointer image);

void MaskImage(vtkSmartPointer<vtkImageData> VTKImage, vtkSmartPointer<vtkImageData> VTKSegmentMask,
               vtkSmartPointer<vtkImageData> VTKMaskedImage);

void WriteWhitePatches(itk::ImageRegion<2> imageRegion, std::vector<itk::ImageRegion<2> > regions, std::string filename);

void VTKImageToITKImage(vtkImageData* image, UnsignedCharScalarImageType::Pointer outputImage);

std::vector<itk::Index<2> > PolyDataToPixelList(vtkPolyData* polydata);

unsigned int CountNonZeroPixels(UnsignedCharScalarImageType::Pointer image);

std::vector<itk::Index<2> > BinaryImageToPixelList(UnsignedCharScalarImageType::Pointer image);

itk::ImageRegion<2> GetCenteredRegionRadius(itk::Index<2> centerPixel, unsigned int radius);

bool IsValidPatch(itk::ImageRegion<2> patch, UnsignedCharScalarImageType::Pointer mask);

bool IsIntersectTargetRegion(itk::ImageRegion<2> patch, UnsignedCharScalarImageType::Pointer mask);

void WritePixels(itk::Size<2> imageSize, std::vector<itk::Index<2> > pixels, std::string filename);

/////////////////// Templated functions (defined in Helpers.txx) /////////////////////////

template <typename TImage>
unsigned int FindClosestIndexGeodesic(typename TImage::Pointer image, std::vector<itk::Index<2> > listOfPixels,
                                      itk::Index<2> queryPixel);

template <typename TImage>
unsigned int CountPixelsWithValue(typename TImage::Pointer image, typename TImage::PixelType pixel);

template <typename TImage>
void ChangeValue(typename TImage::Pointer image, typename TImage::PixelType sourceValue, typename TImage::PixelType targetValue);

template <typename TImage>
std::set<typename TImage::PixelType> GetUniqueValues(typename TImage::Pointer image);

template <typename TImage>
std::set<typename TImage::PixelType> GetNonNegativeUniqueValues(typename TImage::Pointer image);

template <typename TImage>
typename TImage::PixelType GetMaxValue(typename TImage::Pointer image);

template <typename TImage>
std::vector<itk::Index<2> > FindPixelsWithValue(typename TImage::Pointer image, typename TImage::PixelType pixel);

template <typename TImage>
void DeepCopy(typename TImage::Pointer input, typename TImage::Pointer output);

// Convert an ITK image to a VTK image for display
template <typename TImage>
void ITKImageToVTKImage(typename TImage::Pointer image, vtkSmartPointer<vtkImageData> outputImage);

// Specialization for image types with pixel types without [] operator
template <>
void ITKImageToVTKImage<UnsignedCharScalarImageType>(UnsignedCharScalarImageType::Pointer image,
                                                     vtkSmartPointer<vtkImageData> outputImage);

template <typename TPixel>
double SquaredDifference(TPixel, TPixel);

template <typename TImage>
void WriteImage(typename TImage::Pointer image, std::string filename);

template <typename TImage>
void WriteColorMappedImage(typename TImage::Pointer image, std::string filename);

template<typename TImage>
void CastAndWriteImage(typename TImage::Pointer image, std::string filename);

template <typename TImage>
void CopySelfPatchIntoTargetRegion(typename TImage::Pointer image, MaskImageType::Pointer mask,
                                   itk::ImageRegion<2> sourceRegion, itk::ImageRegion<2> targetRegion);

template <typename TImage>
void CopyPatchIntoImage(typename TImage::Pointer sourceImage, typename TImage::Pointer targetImage,
                                   itk::ImageRegion<2> sourceRegion, itk::ImageRegion<2> targetRegion);

template <typename T>
bool ContainsElement(std::vector<T> vec, T element);

template <typename T>
bool HasNeighborWithValue(typename T::Pointer image, itk::Index<2> pixel, typename T::PixelType value);

template <typename T>
unsigned int FindClosestIndexWithMaxGeoDistance(std::vector<itk::Index<2> > listOfPixels, itk::Index<2> queryPixel,
                                                typename T::Pointer image, int maxDistance);

} // end namespace

struct IndexDistance
{
  unsigned int Index;
  float Distance;
};

bool operator<(IndexDistance object1, IndexDistance object2);

#include "Helpers.txx"

#endif
