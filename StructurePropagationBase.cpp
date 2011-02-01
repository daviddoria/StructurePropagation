/*
Copyright (C) 2010 David Doria, daviddoria@gmail.com

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

#include "StructurePropagationBase.h"
#include "Types.h"

#include "itkImageRegionConstIterator.h"
#include "itkRegionOfInterestImageFilter.h"

StructurePropagationBase::StructurePropagationBase()
{
  this->Mask = UnsignedCharScalarImageType::New();

  this->PatchRadius = 5;
}

void StructurePropagationBase::SetMask(UnsignedCharScalarImageType::Pointer mask)
{
  this->Mask->Graft(mask);
}

void StructurePropagationBase::SetPropagationLine(UnsignedCharScalarImageType::Pointer propagationLine)
{
  this->PropagationLine->Graft(propagationLine);
}

void StructurePropagationBase::SetPatchRadius(unsigned int radius)
{
  this->PatchRadius = radius;
}


void StructurePropagationBase::ComputeSourcePatchRegions(std::vector<itk::Index<2> > pixels)
{
  //std::cout << "Patch radius: " << this->PatchRadius << std::endl;

  for(unsigned int i = 0; i < pixels.size(); i++)
    {
    // Create a region centered on the pixels
    itk::Size<2> patchSize;
    patchSize.Fill(this->PatchRadius*2 + 1);

    //std::cout << "Patch size: " << patchSize << std::endl;
    itk::Index<2> corner;
    corner[0] = pixels[i][0] - this->PatchRadius;
    corner[1] = pixels[i][1] - this->PatchRadius;

    itk::ImageRegion<2> region(corner, patchSize);

    // Check if the region intersects the Mask or the edge of the image (both are not allowed)
    itk::ImageRegionConstIterator<UnsignedCharScalarImageType> imageIterator(this->Mask, region);

    bool validRegion = true;
    while(!imageIterator.IsAtEnd())
      {
      if(imageIterator.Get() > 0) // the mask is non-zero
        {
        validRegion = false;
        break;
        }
      ++imageIterator;
      }

    // If it does not, then add the region to the list of source patches
    if(validRegion)
      {
      this->SourcePatchRegions.push_back(region);
      }

    }
}

void StructurePropagationBase::ComputeTargetPatchRegions(std::vector<itk::Index<2> > pixels)
{
  // This function assumes that the stroke starts outside the masked region, enters the masked region, then exits the masked region
  bool enteredMaskedRegion = false;

  for(unsigned int i = 0; i < pixels.size(); i++)
    {
    // Create a region centered on the pixels
    itk::Size<2> patchSize;
    patchSize.Fill(this->PatchRadius*2 + 1);

    itk::Index<2> corner;
    corner[0] = pixels[i][0] - this->PatchRadius;
    corner[1] = pixels[i][1] - this->PatchRadius;

    itk::ImageRegion<2> region(corner, patchSize);

    // Keep patches which intersect the mask by finding the first one, then moving PatchRadius and continuing
    itk::ImageRegionConstIterator<UnsignedCharScalarImageType> imageIterator(this->Mask, region);

    bool inMaskedRegion = false;
    while(!imageIterator.IsAtEnd())
      {
      if(imageIterator.Get() > 0) // the mask is non-zero
        {
        inMaskedRegion = true;
        enteredMaskedRegion = true; // this is a flag that gets set once during the whole loop
        break; // stop checking the current patch
        }
      ++imageIterator;
      }

    // Add the region to the list of target patches
    if(inMaskedRegion)
      {
      this->TargetPatchRegions.push_back(region);
      i+=this->PatchRadius; // skip the next few pixels
      itk::Point<int, 2> point;
      point[0] = pixels[i][0];
      point[1] = pixels[i][1];
      //NodeLocations.push_back(pixels[i]);
      NodeLocations.push_back(point);
      }

    // Stop when we find the first patch that is not in the target region after we have gone through the hole
    if(!inMaskedRegion && enteredMaskedRegion)
      {
      break; // we are done
      }

    }
}


std::vector<itk::ImageRegion<2> > StructurePropagationBase::GetSourcePatchRegions()
{
  return this->SourcePatchRegions;
}

std::vector<itk::ImageRegion<2> > StructurePropagationBase::GetTargetPatchRegions()
{
  return this->TargetPatchRegions;
}

void StructurePropagationBase::ExtractRegionsOfPropagationPath()
{
  // Extract the regions of the propagation path
  this->SourceStrokePaths.clear();
  typedef itk::RegionOfInterestImageFilter< UnsignedCharScalarImageType, UnsignedCharScalarImageType > ROIFilterType;

  for(unsigned int i = 0; i < this->SourcePatchRegions.size(); i++)
    {
    ROIFilterType::Pointer roiFilter = ROIFilterType::New();
    roiFilter->SetRegionOfInterest(this->SourcePatchRegions[i]);
    roiFilter->SetInput(this->PropagationLine);
    this->SourceStrokePaths.push_back(roiFilter->GetOutput());
    }

  this->TargetStrokePaths.clear();
  for(unsigned int i = 0; i < this->TargetPatchRegions.size(); i++)
    {
    ROIFilterType::Pointer roiFilter = ROIFilterType::New();
    roiFilter->SetRegionOfInterest(this->TargetPatchRegions[i]);
    roiFilter->SetInput(this->PropagationLine);
    this->TargetStrokePaths.push_back(roiFilter->GetOutput());
    }
}

void StructurePropagationBase::ComputePatchRegions(std::vector<itk::Index<2> > pixels)
{
  this->ComputeSourcePatchRegions(pixels);
  this->ComputeTargetPatchRegions(pixels);
}