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
void StructurePropagation<TImage>::SetSourceRegions(const std::vector<itk::ImageRegion<2> >& sourceRegions)
{
  this->SourceRegions = sourceRegions;
}

template <typename TImage>
void StructurePropagation<TImage>::SetTargetRegions(const std::vector<itk::ImageRegion<2> >& targetRegions)
{
  this->TargetRegions = targetRegions;
}

template <typename TImage>
void StructurePropagation<TImage>::PropagateStructure()
{
  std::cout << "StructurePropagation: PropagateStructure()" << std::endl;
  std::cout << "There are " << this->SourceRegions.size() << " source patch regions (labels)." << std::endl;
  std::cout << "There are " << this->TargetRegions.size() << " target patch regions (nodes)." << std::endl;

  StructurePropagationDynamicProgramming<TImage> dynamicProgramming;
  dynamicProgramming.SetImage(this->Image);
  dynamicProgramming.SetMask(this->MaskImage);
  dynamicProgramming.SetLabelSet(this->SourceRegions);
  dynamicProgramming.SetNodes(this->TargetRegions);
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
                                                        this->SourceRegions[solution[i]],
                                                        this->TargetRegions[i]);
    tempMask->SetValid(this->TargetRegions[i]);
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
void StructurePropagation<TImage>::ComputeSourceRegions(const PropagationLineImageType* const propagationLineOutsideMask)
{
  // Start fresh
  this->SourceRegions.clear();

  std::vector<itk::Index<2> > nonZeroPixels = ITKHelpers::GetNonZeroPixels(propagationLineOutsideMask);

  for(unsigned int i = 0; i < nonZeroPixels.size(); ++i)
  {
    itk::ImageRegion<2> region = ITKHelpers::GetRegionInRadiusAroundPixel(nonZeroPixels[i], this->PatchRadius);
    if(this->MaskImage->GetLargestPossibleRegion().IsInside(region))
    {
      if(this->MaskImage->IsValid(region))
      {
        this->SourceRegions.push_back(region);
      }
    }
  }
}

template <typename TImage>
void StructurePropagation<TImage>::ComputeTargetRegionsOpenContour(const PropagationLineImageType* const propagationLineInsideMask)
{
  // Start fresh
  this->TargetRegions.clear();

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
    this->TargetRegions.push_back(region);
  }
}


template <typename TImage>
void StructurePropagation<TImage>::ComputeTargetRegionsClosedContour(const PropagationLineImageType* const propagationLineInsideMask)
{
  ITKHelpers::WriteImage(propagationLineInsideMask, "propagation.png");

  // Start fresh
  this->TargetRegions.clear();

  itk::Index<2> pixelOnContour = ITKHelpers::FindFirstNonZeroPixel(propagationLineInsideMask);
  std::vector<itk::Index<2> > orderedPixels = ITKHelpers::GetClosedContourOrdering(propagationLineInsideMask, pixelOnContour);

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
    this->TargetRegions.push_back(region);
  }
}

template <typename TImage>
std::vector<itk::ImageRegion<2> > StructurePropagation<TImage>::GetSourcePatchRegions()
{
  return this->SourceRegions;
}

template <typename TImage>
std::vector<itk::ImageRegion<2> > StructurePropagation<TImage>::GetTargetPatchRegions()
{
  return this->TargetRegions;
}

template <typename TImage>
void StructurePropagation<TImage>::ComputeDonutRegionsPropagationAroundHole()
{
  Mask::BoundaryImageType::Pointer boundaryImage = Mask::BoundaryImageType::New();
  this->MaskImage->FindBoundary(boundaryImage, Mask::VALID);
  ComputeTargetRegionsClosedContour(boundaryImage);

  Mask::Pointer expandedMask = Mask::New();
  expandedMask->DeepCopyFrom(this->MaskImage);
  unsigned int expandFactor = 1;
  expandedMask->ExpandHole(this->PatchRadius * expandFactor);

  typedef itk::XorImageFilter<Mask> XORFilterType;
  XORFilterType::Pointer xorFilter = XORFilterType::New();
  xorFilter->SetInput1(this->MaskImage);
  xorFilter->SetInput2(expandedMask);
  xorFilter->Update();

  Mask::Pointer donutMask = Mask::New();
  donutMask->DeepCopyFrom(xorFilter->GetOutput());
  donutMask->CopyInformationFrom(this->MaskImage);

  ITKHelpersTypes::UnsignedCharScalarImageType::Pointer donutImage = ITKHelpersTypes::UnsignedCharScalarImageType::New();
  ITKHelpers::InitializeImage(donutImage.GetPointer(), this->Image->GetLargestPossibleRegion());
  ITKHelpers::SetRegionToConstant(donutImage.GetPointer(), donutImage->GetLargestPossibleRegion(), 255);
  donutMask->ApplyToScalarImage(donutImage.GetPointer(), 0);

  ITKHelpers::WriteImage(donutImage.GetPointer(), "donut.png");

  this->ComputeSourceRegions(donutImage);
}


template <typename TImage>
void StructurePropagation<TImage>::ComputeRegionsPropagationAroundHole()
{
  // Computing a donut of source regions produces way too many source regions. Instead, just keep the outer ring of pixels
  // around the donut.
  Mask::BoundaryImageType::Pointer boundaryImage = Mask::BoundaryImageType::New();
  this->MaskImage->FindBoundary(boundaryImage, Mask::VALID);
  ComputeTargetRegionsClosedContour(boundaryImage);

  Mask::Pointer expandedMask = Mask::New();
  expandedMask->DeepCopyFrom(this->MaskImage);
  unsigned int expandFactor = 1;
  expandedMask->ExpandHole(this->PatchRadius * expandFactor);

  Mask::BoundaryImageType::Pointer expandedBoundaryImage = Mask::BoundaryImageType::New();
  expandedMask->FindBoundary(expandedBoundaryImage, Mask::VALID);

  ITKHelpers::WriteImage(expandedBoundaryImage.GetPointer(), "expandedBoundaryImage.png");

  this->ComputeSourceRegions(expandedBoundaryImage);
}

template <typename TImage>
void StructurePropagation<TImage>::ComputeRegionsPropagationThroughHole()
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

  this->ComputeSourceRegions(propagationLineOutsideMask);
  this->ComputeTargetRegionsOpenContour(propagationLineInsideMask);
}

template <typename TImage>
void StructurePropagation<TImage>::SetPropagationLineImage(PropagationLineImageType*
                                                           const propagationLineImage)
{
  ITKHelpers::DeepCopy(propagationLineImage, this->PropagationLineImage.GetPointer());
}

#endif
