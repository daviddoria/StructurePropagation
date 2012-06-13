#ifndef StructurePropagationDynamicProgramming_hpp
#define StructurePropagationDynamicProgramming_hpp

#include "StructurePropagationDynamicProgramming.h"

#include "Mask/ITKHelpers/ITKHelpers.h"

void StructurePropagationDynamicProgramming::SetImage(ImageType* const image)
{
  ITKHelpers::DeepCopy(image, this->Image.GetPointer());
}

float StructurePropagationDynamicProgramming::BinaryEnergy(const LabelType& labelA, const unsigned int nodeA,
                                           const LabelType& labelB, const unsigned int nodeB)
{
  itk::ImageRegion<2> targetRegion1 = this->TargetPatchRegions[node1];
  itk::ImageRegion<2> targetRegion2 = this->TargetPatchRegions[node2];

  itk::ImageRegion<2> sourceRegion1 = this->SourcePatchRegions[label1];
  itk::ImageRegion<2> sourceRegion2 = this->SourcePatchRegions[label2];

  // Find the overlapping region
  itk::ImageRegion<2> overlap = targetRegion1;
  overlap.Crop(targetRegion2);

  // Find the offset of the corner of the overlap region from each of the patches
  itk::Offset<2> offset1 = overlap.GetIndex() - targetRegion1.GetIndex();
  itk::Offset<2> offset2 = overlap.GetIndex() - targetRegion2.GetIndex();

  itk::ImageRegion<2> overlapRegion1(sourceRegion1.GetIndex() + offset1, overlap.GetSize());
  itk::ImageRegion<2> overlapRegion2(sourceRegion2.GetIndex() + offset2, overlap.GetSize());

  // sum of the normalized squared differences in the region where the two patches overlap
  itk::ImageRegionConstIterator<TImage> patch1Iterator(this->Image, overlapRegion1);
  itk::ImageRegionConstIterator<TImage> patch2Iterator(this->Image, overlapRegion2);

  double ssd = 0; // sum of squared differences
  while(!patch1Iterator.IsAtEnd())
    {
    // Get the value of the current pixel
    typename TImage::PixelType pixel1 = patch1Iterator.Get();
    typename TImage::PixelType pixel2 = patch2Iterator.Get();
    ssd += Helpers::SquaredDifference(pixel1, pixel2);

    ++patch1Iterator;
    ++patch2Iterator;
    }

  /*
  std::cout << "Binary cost for nodes " << node1 << " and " << node2
            << " with labels " << label1 << " and " << label2 << " is " << ssd << std::endl;
  */

  if(ssd <= 0)
    {
    return 0.0;
    }
  else
    {
    return ssd/static_cast<double>(overlapRegion1.GetNumberOfPixels());
    }
}

float StructurePropagationDynamicProgramming::UnaryEnergy(const LabelType& label, const unsigned int node)
{
  // "sum of the normalized squared differences in the region of the patch which overlaps the source region"
  itk::ImageRegionConstIterator<TImage> imageIterator(this->Image, this->TargetPatchRegions[node]);
  itk::ImageRegionConstIterator<UnsignedCharScalarImageType> maskIterator(this->Mask,this->TargetPatchRegions[node]);

  itk::ImageRegionConstIterator<TImage> patchIterator(this->Image, this->SourcePatchRegions[label]);

  double ssd = 0; // sum of squared differences
  unsigned int numberOfPixels = 0;
  while(!imageIterator.IsAtEnd())
    {
    if(maskIterator.Get() == 0) // this is the region we want (assuming mask is black in the source region)
      {
      // Get the value of the current pixel
      typename TImage::PixelType pixel1 = imageIterator.Get();
      typename TImage::PixelType pixel2 = patchIterator.Get();
      ssd += Helpers::SquaredDifference(pixel1, pixel2);
      numberOfPixels++;
      }

    ++imageIterator;
    ++maskIterator;
    ++patchIterator;
    }
  if(numberOfPixels <= 0)
    {
    return 0;
    }
  else
    {
    return ssd/static_cast<double>(numberOfPixels);
    }
}

#endif
