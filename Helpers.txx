#include "itkCastImageFilter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkScalarToRGBColormapImageFilter.h"

namespace Helpers
{

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
  rgbfilter->SetColormap( RGBFilterType::Hot );
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

} // end namespace