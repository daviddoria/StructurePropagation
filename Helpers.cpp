#include "Helpers.h"
#include "Types.h"

#include "itkImageRegionIterator.h"
#include "itkBresenhamLine.h"
#include "itkImageRegionConstIterator.h"

#include <vtkPolyData.h>

void IndicesToBinaryImage(std::vector<itk::Index<2> > indices, UnsignedCharScalarImageType::Pointer image)
{
  //std::cout << "Setting " << indices.size() << " points to non-zero." << std::endl;

  // Blank the image
  itk::ImageRegionIterator<UnsignedCharScalarImageType> imageIterator(image,image->GetLargestPossibleRegion());
  while(!imageIterator.IsAtEnd())
    {
    imageIterator.Set(0);
    ++imageIterator;
    }

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
void ITKImagetoVTKImage<UnsignedCharScalarImageType>(UnsignedCharScalarImageType::Pointer image, vtkImageData* outputImage)
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