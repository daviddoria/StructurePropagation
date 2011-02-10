#ifndef HELPERS_H
#define HELPERS_H

#include "Types.h"

#include "itkImage.h"
#include "itkIndex.h"

#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkPolyData.h>
#include <vtkUnsignedCharArray.h>

#include <string>

namespace Helpers
{

void DrawLineInImage(UnsignedCharScalarImageType::Pointer image, itk::Index<2> p0, itk::Index<2> p1);

void NumberPixels(UnsignedCharScalarImageType::Pointer binaryImage, IntScalarImageType::Pointer numberedPixels);

itk::Index<2> GetCorner(itk::Index<2> center, unsigned int radius);

void GetDilatedImage(UnsignedCharScalarImageType::Pointer inputImage, UnsignedCharScalarImageType::Pointer dilatedImage);

double SumOfShortestDistances(std::vector<itk::Index<2> > indices1, std::vector<itk::Index<2> > indices2);
double SumOfShortestDistances(std::vector<itk::Point<int, 2> > points1, std::vector<itk::Point<int, 2> > points2);

// Mark each pixel at the specified 'indices' as a non-zero pixel in 'image'
void IndicesToBinaryImage(std::vector<itk::Index<2> > indices, UnsignedCharScalarImageType::Pointer image);

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

// Convert an ITK image to a VTK image for display
template <typename TImageType>
void ITKImageToVTKImage(typename TImageType::Pointer image, vtkSmartPointer<vtkImageData> outputImage);

// Specialization for image types with pixel types without [] operator
template <>
void ITKImageToVTKImage<UnsignedCharScalarImageType>(UnsignedCharScalarImageType::Pointer image,
                                                     vtkSmartPointer<vtkImageData> outputImage);

template <typename TPixelType>
double SquaredDifference(TPixelType, TPixelType);

template <typename TImageType>
void WriteImage(typename TImageType::Pointer image, std::string filename);

template <typename TImageType>
void WriteColorMappedImage(typename TImageType::Pointer image, std::string filename);

template<typename TImageType>
void CastAndWriteImage(typename TImageType::Pointer image, std::string filename);

template <typename TImageType>
void CopySelfPatchIntoTargetRegion(typename TImageType::Pointer image, MaskImageType::Pointer mask, itk::ImageRegion<2> sourceRegion, itk::ImageRegion<2> targetRegion);

} // end namespace

#include "Helpers.txx"

#endif
