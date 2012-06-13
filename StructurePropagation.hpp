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

#ifndef StructurePropagation_HPP
#define StructurePropagation_HPP

// Appease syntax parser
#include "StructurePropagation.h"

// ITK
#include "itkBinaryNotImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkContourMeanDistanceImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIterator.h"
#include "itkMaskImageFilter.h"
#include "itkPasteImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkXorImageFilter.h"

// Submodules
#include "DynamicProgramming/Helpers/Helpers.h"
#include "Mask/ITKHelpers/ITKHelpers.h"
#include "Mask/MaskOperations.h"

// Custom
#include "StructurePropagationDynamicProgramming.h"

template <typename TImage>
StructurePropagation<TImage>::StructurePropagation()
{
  // Initializations
  this->Image = TImage::New();
  this->OutputImage = TImage::New();
  this->MaskImage = Mask::New();
  this->PropagationLineImage = PropagationLineImageType::New();

  // Default Values
  this->PatchRadius = 10;
}

template <typename TImage>
TImage* StructurePropagation<TImage>::GetOutputImage()
{
  return this->OutputImage.GetPointer();
}

template <typename TImage>
void StructurePropagation<TImage>::SetImage(TImage* const image)
{
  ITKHelpers::DeepCopy(image, this->Image.GetPointer());
}

template <typename TImage>
void StructurePropagation<TImage>::ClearEverything()
{
  this->SourcePatchRegions.clear();
  this->TargetPatchRegions.clear();
}

template <typename TImage>
void StructurePropagation<TImage>::PropagateStructure()
{
  ComputePatchRegions();

  std::cout << "There are " << this->SourcePatchRegions.size() << " source patch regions (labels)." << std::endl;
  std::cout << "There are " << this->TargetPatchRegions.size() << " target patch regions (nodes)." << std::endl;

  StructurePropagationDynamicProgramming<TImage> dynamicProgramming;
  dynamicProgramming.SetImage(this->Image);
  dynamicProgramming.SetMask(this->MaskImage);
  dynamicProgramming.SetLabelSet(this->SourcePatchRegions);
  dynamicProgramming.SetNodes(this->TargetPatchRegions);
  std::vector<unsigned int> solution = dynamicProgramming.Optimize();

  {// Debug only
  std::cout << "Solution labels: ";
  for(unsigned int i = 0; i < solution.size(); i++)
    {
    std::cout << solution[i] << " ";
    }
  std::cout << std::endl;
  }

  ITKHelpers::DeepCopy(this->Image.GetPointer(), this->OutputImage.GetPointer());

  Mask::Pointer tempMask = Mask::New();
  tempMask->DeepCopyFrom(this->MaskImage);

  // Fill in the image with the best patches
  for(unsigned int i = 0; i < solution.size(); i++)
    {
    MaskOperations::CopySelfPatchIntoHoleOfTargetRegion(this->OutputImage.GetPointer(), tempMask,
                                                        this->SourcePatchRegions[solution[i]],
                                                        this->TargetPatchRegions[i]);
    tempMask->SetValid(this->TargetPatchRegions[i]);
    }

  std::cout << "Finished propagating structure!" << std::endl;
}


template <typename TImage>
void StructurePropagation<TImage>::SetMask(Mask* const mask)
{
  this->MaskImage->DeepCopyFrom(mask);
}

template <typename TImage>
void StructurePropagation<TImage>::SetPatchRadius(unsigned int radius)
{
  this->PatchRadius = radius;
}

template <typename TImage>
void StructurePropagation<TImage>::ComputeSourcePatchRegions(const PropagationLineImageType* const propagationLineOutsideMask)
{
  // Start fresh
  this->SourcePatchRegions.clear();

  std::vector<itk::Index<2> > nonZeroPixels = ITKHelpers::GetNonZeroPixels(propagationLineOutsideMask);

  for(unsigned int i = 0; i < nonZeroPixels.size(); ++i)
  {
    itk::ImageRegion<2> region = ITKHelpers::GetRegionInRadiusAroundPixel(nonZeroPixels[i], this->PatchRadius);
    if(this->MaskImage->IsValid(region))
    {
      this->SourcePatchRegions.push_back(region);
    }
  }
}

template <typename TImage>
void StructurePropagation<TImage>::ComputeTargetPatchRegions(const PropagationLineImageType* const propagationLineInsideMask)
{
  // Start fresh
  this->TargetPatchRegions.clear();

  itk::Index<2> pixelOnContour = ITKHelpers::FindFirstNonZeroPixel(propagationLineInsideMask);
  std::vector<itk::Index<2> > orderedPixels = ITKHelpers::GetOpenContourOrdering(propagationLineInsideMask, pixelOnContour);

  // Sample the pixels on the contour so they are about a patchRadius away from each other
  std::vector<itk::Index<2> > sampledPixels;

  // The first pixel is definitely in the list (should be on the hole boundary)
  sampledPixels.push_back(orderedPixels[0]);

  for(unsigned int i = 1; i < orderedPixels.size(); i++) // start at 1 because we have already used 0
  {
    if(ITKHelpers::IndexDistance(sampledPixels[sampledPixels.size() - 1], orderedPixels[i]) >= this->PatchRadius)
    {
      sampledPixels.push_back(orderedPixels[i]);
    }
  }

  // Always add the last pixel (if it wasn't already used)
  if(sampledPixels[sampledPixels.size() - 1] != orderedPixels[orderedPixels.size() - 1])
  {
    sampledPixels.push_back(orderedPixels[orderedPixels.size() - 1]);
  }

  { // debug only
  ITKHelpersTypes::UnsignedCharScalarImageType::Pointer centersImage = ITKHelpersTypes::UnsignedCharScalarImageType::New();
  ITKHelpers::InitializeImage(centersImage.GetPointer(), this->Image->GetLargestPossibleRegion());
  ITKHelpers::SetPixels(centersImage.GetPointer(), sampledPixels, 255);
  ITKHelpers::WriteImage(centersImage.GetPointer(), "nodeCenters.png");
  }
  // Create regions centered on the pixels
  for(unsigned int i = 0; i < sampledPixels.size(); ++i)
  {
    itk::ImageRegion<2> region = ITKHelpers::GetRegionInRadiusAroundPixel(sampledPixels[i], this->PatchRadius);
    this->TargetPatchRegions.push_back(region);
  }
}

template <typename TImage>
std::vector<itk::ImageRegion<2> > StructurePropagation<TImage>::GetSourcePatchRegions()
{
  return this->SourcePatchRegions;
}

template <typename TImage>
std::vector<itk::ImageRegion<2> > StructurePropagation<TImage>::GetTargetPatchRegions()
{
  return this->TargetPatchRegions;
}

template <typename TImage>
void StructurePropagation<TImage>::ComputePatchRegions()
{
  // Extract the part of the lines in the source region
  PropagationLineImageType::Pointer propagationLineOutsideMask = PropagationLineImageType::New();
  ITKHelpers::DeepCopy(this->PropagationLineImage.GetPointer(), propagationLineOutsideMask.GetPointer());
  this->MaskImage->ApplyToScalarImage(propagationLineOutsideMask.GetPointer(), 0);

  ITKHelpers::WriteImage(propagationLineOutsideMask.GetPointer(), "propagationLineOutsideMask.png");

  // Extract the part of the lines in the target region
  typedef itk::XorImageFilter<PropagationLineImageType> XORFilterType;
  XORFilterType::Pointer xorFilter = XORFilterType::New();
  xorFilter->SetInput1(propagationLineOutsideMask);
  xorFilter->SetInput2(this->PropagationLineImage);
  xorFilter->Update();

  PropagationLineImageType::Pointer propagationLineInsideMask = PropagationLineImageType::New();
  ITKHelpers::DeepCopy(xorFilter->GetOutput(), propagationLineInsideMask.GetPointer());

  ITKHelpers::WriteImage(propagationLineInsideMask.GetPointer(), "propagationLineInsideMask.png");

  this->ComputeSourcePatchRegions(propagationLineOutsideMask);
  this->ComputeTargetPatchRegions(propagationLineInsideMask);
}

template <typename TImage>
void StructurePropagation<TImage>::SetPropagationLineImage(PropagationLineImageType*
                                                           const propagationLineImage)
{
  ITKHelpers::DeepCopy(propagationLineImage, this->PropagationLineImage.GetPointer());
}

#endif
