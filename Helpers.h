#ifndef HELPERS_H
#define HELPERS_H

#include "itkImage.h"
#include "itkIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"

#include "vtkSmartPointer.h"
#include "vtkImageData.h"

#include "Types.h"

class vtkPolyData;

// Mark each pixel at the specified 'indices' as a non-zero pixel in 'image'
void IndicesToBinaryImage(std::vector<itk::Index<2> > indices, UnsignedCharScalarImageType::Pointer image);

void MaskImage(vtkSmartPointer<vtkImageData> VTKImage, vtkSmartPointer<vtkImageData> VTKSegmentMask, vtkSmartPointer<vtkImageData> VTKMaskedImage);

// Convert an ITK image to a VTK image for display
template <typename TImageType>
void ITKImagetoVTKImage(typename TImageType::Pointer image, vtkImageData* outputImage)
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

// Specialization for image types with pixel types without [] operator
template <>
void ITKImagetoVTKImage<UnsignedCharScalarImageType>(UnsignedCharScalarImageType::Pointer image, vtkImageData* outputImage);

std::vector<itk::Index<2> > PolyDataToPixelList(vtkPolyData* polydata);

std::vector<itk::Index<2> > BinaryImageToPixelList(UnsignedCharScalarImageType::Pointer image);

#endif