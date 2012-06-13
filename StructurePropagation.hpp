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

template <typename TImage>
StructurePropagation<TImage>::StructurePropagation()
{
  // Initializations
  this->Image = TImage::New();
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
  this->NodeLocations.clear();
  this->Edges.clear();
}

template <typename TImage>
void StructurePropagation<TImage>::PropagateStructure()
{
  // Find intersections of target line with hole boundary and add them as
  // additional intersections (so nodes are required to be placed there)
  for(IndexSetType::iterator iterator = this->PropagationLine.begin();
      iterator != this->PropagationLine.end(); iterator++)
    {
    if(this->MaskImage->IsHole(*iterator)) // the pixel is in the target region
      {
      if(ITKHelpers::HasNeighborWithValue(this->MaskImage, *iterator, 0))
        {
        this->PropagationLineIntersections.insert(*iterator);
        std::cout << "Inserted boundary intersection." << std::endl;
        }
      }
    }

  // Create the graphical model
  unsigned int numberOfNodes = this->TargetPatchRegions.size();
  unsigned int numberOfLabels = this->SourcePatchRegions.size();
  std::cout << "PropagateStructure: There are " << numberOfNodes << " nodes." << std::endl;
  std::cout << "PropagateStructure: There are " << numberOfLabels << " labels." << std::endl;

  // Fill in the image with the best patches
//   for(unsigned int i = 0; i < numberOfNodes; i++)
//     {
//     itk::ImageRegion<2> destinationRegion =
//            ITKHelpers::GetRegionInRadiusAroundPixel(this->NodeLocations[i], this->PatchRadius);
//     MaskOperations::CopySelfPatchIntoHoleOfTargetRegion(this->OutputImage.GetPointer(), this->MaskImage,
//                                                         this->SourcePatchRegions[result[i]],
//                                                         destinationRegion);
//     }

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

//   ITKHelpersTypes::UnsignedCharScalarImageType::Pointer dilatedLineImage =
//            ITKHelpersTypes::UnsignedCharScalarImageType::New();
//   ITKHelpers::DilateImage(this->PropagationLineSourceImage, this->SourceLineWidth/2, dilatedLineImage);
//   Helpers::WriteImage(dilatedLineImage, "DilatedSourceLines.png");

  // The dilated source line image will bleed into the target region,
  // but this is ok because we ensure that only patches completely in the sourcce region are kept

}

template <typename TImage>
void StructurePropagation<TImage>::ComputeTargetPatchRegions(const PropagationLineImageType* const propagationLineInsideMask)
{
  // Start fresh
  this->TargetPatchRegions.clear();

  itk::Size<2> patchSize;
  patchSize.Fill(this->PatchRadius*2 + 1);

  for(unsigned int i = 0; i < this->NodeLocations.size(); i++)
    {
    // Create a region centered on the pixels
    itk::Index<2> region = ITKHelpers::GetRegionInRadiusAroundPixel(this->NodeLocations[i], this->PatchRadius);

    // Keep patches which intersect the target region
//     if(Helpers::IsIntersectTargetRegion(region, this->Mask));
//       {
//       this->TargetPatchRegions.push_back(region);
//       }
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
  this->MaskImage->ApplyToImage(propagationLineOutsideMask, 0);

  ITKHelpers::WriteImage(propagationLineOutsideMask.GetPointer(), "propagationLineOutsideMask.png");

  // Extract the part of the lines in the target region
  PropagationLineImageType::Pointer propagationLineInsideMask = PropagationLineImageType::New();
  ITKHelpers::DeepCopy(this->PropagationLineImage.GetPointer(), propagationLineInsideMask.GetPointer());
  typedef itk::XorImageFilter<PropagationLineImageType> XORFilterType;
  XORFilterType::Pointer xorFilter = XORFilterType::New();
  xorFilter->SetInput1(propagationLineOutsideMask);
  xorFilter->SetInput1(this->PropagationLineImage);
  xorFilter->Update();

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
