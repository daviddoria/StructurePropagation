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
template <typename TImage>
unsigned int CountPixelsWithValue(typename TImage::Pointer image, typename TImage::PixelType pixel);

unsigned int CountNumberOfNonNegativeValues(IntScalarImageType::Pointer image);

template <typename TImage>
typename TImage::PixelType GetMaxValue(typename TImage::Pointer image);

template <typename TImage>
std::vector<itk::Index<2> > FindPixelsWithValue(typename TImage::Pointer image, typename TImage::PixelType pixel);

std::vector<itk::Index<2> > GetLabelCenters(IntScalarImageType::Pointer image);
itk::Index<2> AverageIndex(std::vector<itk::Index<2> > pixels);

void DrawLineInImage(UnsignedCharScalarImageType::Pointer image, itk::Index<2> p0, itk::Index<2> p1);

void NumberPixels(UnsignedCharScalarImageType::Pointer binaryImage, IntScalarImageType::Pointer numberedPixels);
void ClusterNumberedPixels(IntScalarImageType::Pointer numberedPixels, IntScalarImageType::Pointer clusteredNumberedPixels);
void RelabelSequential(IntScalarImageType::Pointer numberedPixels);

template <typename TImage>
void DeepCopy(typename TImage::Pointer input, typename TImage::Pointer output);

itk::Index<2> GetCorner(itk::Index<2> center, unsigned int radius);

std::vector<itk::Offset<2> > CreateOffsetsInRadius(itk::Size<2> radius);

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

void WritePixels(itk::Size<2> imageSize, std::vector<itk::Index<2> > pixels, std::string filename);

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

} // end namespace

#include "Helpers.txx"

#endif
