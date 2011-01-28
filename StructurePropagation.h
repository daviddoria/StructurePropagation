#ifndef StructurePropagation_H
#define StructurePropagation_H

#include "StructurePropagationBase.h"
#include "itkImageRegionConstIterator.h"

#include "Types.h"

template <typename TImageType>
class StructurePropagation : public StructurePropagationBase
{
public:
  StructurePropagation(typename TImageType::Pointer image);

  void SetImage(typename TImageType::Pointer image);

  void ExtractSourcePatches(std::vector<itk::Index<2> > pixels);
  void ComputeTargetPatches(std::vector<itk::Index<2> > pixels);

private:
  typename TImageType::Pointer Image;

};

template <typename TImageType>
StructurePropagation<TImageType>::StructurePropagation(typename TImageType::Pointer image)
{
  this->Image = TImageType::New();
  this->Image->Graft(image);
}

template <typename TImageType>
void StructurePropagation<TImageType>::SetImage(typename TImageType::Pointer image)
{
  //this->Image = image;
  this->Image->Graft(image);
}

template <typename TImageType>
void StructurePropagation<TImageType>::ExtractSourcePatches(std::vector<itk::Index<2> > pixels)
{
  for(unsigned int i = 0; i < pixels.size(); i++)
    {
    // Create a region centered on the pixels
    itk::Size<2> patchSize;
    patchSize.Fill(this->PatchRadius*2 + 1);

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
      this->PatchRegions.push_back(region);
      }

    }
}


template <typename TImageType>
void StructurePropagation<TImageType>::ComputeTargetPatches(std::vector<itk::Index<2> > pixels)
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

    // If it does not, then add the region to the list of source patches
    if(inMaskedRegion)
      {
      this->PatchRegions.push_back(region);
      i+=this->PatchRadius; // skip the next few pixels
      }

    if(!inMaskedRegion && enteredMaskedRegion)
      {
      break; // we are done
      }

    }
}

#endif