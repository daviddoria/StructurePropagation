#ifndef HELPERS_H
#define HELPERS_H

#include "itkImage.h"
#include "itkIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"

#include "vtkSmartPointer.h"
#include "vtkImageData.h"

#include "Types.h"

class vtkPolyData;

namespace Helpers
{

void GetDilatedImage(UnsignedCharScalarImageType::Pointer inputImage, UnsignedCharScalarImageType::Pointer dilatedImage);

template <typename TPixelType>
double difference(TPixelType, TPixelType);

double SumOfShortestDistances(std::vector<itk::Index<2> > indices1, std::vector<itk::Index<2> > indices2);
double SumOfShortestDistances(std::vector<itk::Point<int, 2> > points1, std::vector<itk::Point<int, 2> > points2);

// Mark each pixel at the specified 'indices' as a non-zero pixel in 'image'
void IndicesToBinaryImage(std::vector<itk::Index<2> > indices, UnsignedCharScalarImageType::Pointer image);

void MaskImage(vtkSmartPointer<vtkImageData> VTKImage, vtkSmartPointer<vtkImageData> VTKSegmentMask, vtkSmartPointer<vtkImageData> VTKMaskedImage);

void WriteWhitePatches(itk::ImageRegion<2> imageRegion, std::vector<itk::ImageRegion<2> > regions, std::string filename);

void VTKImageToITKImage(vtkImageData* image, UnsignedCharScalarImageType::Pointer outputImage);

// Convert an ITK image to a VTK image for display
template <typename TImageType>
void ITKImageToVTKImage(typename TImageType::Pointer image, vtkSmartPointer<vtkImageData> outputImage);

// Specialization for image types with pixel types without [] operator
template <>
void ITKImageToVTKImage<UnsignedCharScalarImageType>(UnsignedCharScalarImageType::Pointer image, vtkSmartPointer<vtkImageData> outputImage);

std::vector<itk::Index<2> > PolyDataToPixelList(vtkPolyData* polydata);

unsigned int CountNonZeroPixels(UnsignedCharScalarImageType::Pointer image);

std::vector<itk::Index<2> > BinaryImageToPixelList(UnsignedCharScalarImageType::Pointer image);

#include "Helpers.txx"

} // end namespace

#endif
