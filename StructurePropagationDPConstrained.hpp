#ifndef StructurePropagationDPConstrained_hpp
#define StructurePropagationDPConstrained_hpp

// Appease syntax parser
#include "StructurePropagationDPConstrained.h"

// Submodules
#include "Mask/ITKHelpers/ITKHelpers.h"

// ITK
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"

template <typename TImage>
StructurePropagationDPConstrained<TImage>::StructurePropagationDPConstrained()
{
  this->Image = NULL;
  this->MaskImage = NULL;
}

template <typename TImage>
void StructurePropagationDPConstrained<TImage>::SetUnaryBinaryWeight(const float unaryBinaryWeight)
{
  this->UnaryBinaryWeight = unaryBinaryWeight;
}

template <typename TImage>
void StructurePropagationDPConstrained<TImage>::SetMaxNodeLabelDistance(const unsigned int maxNodeLabelDistance)
{
  this->MaxNodeLabelDistance = maxNodeLabelDistance;
}

template <typename TImage>
void StructurePropagationDPConstrained<TImage>::SetImage(TImage* const image)
{
  //ITKHelpers::DeepCopy(image, this->Image.GetPointer());
  this->Image = image;
}

template <typename TImage>
void StructurePropagationDPConstrained<TImage>::SetMask(Mask* const mask)
{
  //this->MaskImage->DeepCopyFrom(mask);
  this->MaskImage = mask;
}

template <typename TImage>
float StructurePropagationDPConstrained<TImage>::BinaryEnergy(const LabelType& labelA, const NodeType& nodeA,
                                                              const LabelType& labelB, const NodeType& nodeB)
{
  // This function assumes that both label regions are entirely valid
//   assert(this->Image->GetLargestPossibleRegion().GetSize()[0] != 0);
//   assert(this->MaskImage->GetLargestPossibleRegion().GetSize()[0] != 0);
  assert(this->Image);
  assert(this->MaskImage);
  assert(this->MaskImage->IsValid(labelA));
  assert(this->MaskImage->IsValid(labelB));

  if(ITKHelpers::IndexDistance(ITKHelpers::GetRegionCenter(labelA), ITKHelpers::GetRegionCenter(nodeA)) > this->MaxNodeLabelDistance ||
     ITKHelpers::IndexDistance(ITKHelpers::GetRegionCenter(labelB), ITKHelpers::GetRegionCenter(nodeB)) > this->MaxNodeLabelDistance)
  {
    return std::numeric_limits<float>::infinity();
  }

//   if(ITKHelpers::IndexDistance(ITKHelpers::GetRegionCenter(labelA), ITKHelpers::GetRegionCenter(nodeA)) > 2 * labelA.GetSize()[0] ||
//      ITKHelpers::IndexDistance(ITKHelpers::GetRegionCenter(labelB), ITKHelpers::GetRegionCenter(nodeB)) > 2 * labelB.GetSize()[0]
//   )
//   {
//     return std::numeric_limits<float>::infinity();
//   }

  // Find the overlapping region
  itk::ImageRegion<2> overlap = nodeA;
  overlap.Crop(nodeB);

  if(overlap.GetNumberOfPixels() == 0)
  {
    std::stringstream ss;
    ss << "There are 0 pixels in the overlap of nodes " << nodeA << " and " << nodeB;
    std::runtime_error(ss.str());
    return 0.0f;
  }

  // Find the offset of the corner of the overlap region from each of the node regions
  itk::Offset<2> nodeA_offset = overlap.GetIndex() - nodeA.GetIndex();
  itk::Offset<2> nodeB_offset = overlap.GetIndex() - nodeB.GetIndex();

  // Get the corresponding regions in the label regions
  itk::ImageRegion<2> labelA_OverlapRegion(labelA.GetIndex() + nodeA_offset, overlap.GetSize());
  itk::ImageRegion<2> labelB_OverlapRegion(labelB.GetIndex() + nodeB_offset, overlap.GetSize());

  // sum of the squared differences in the region where the two patches overlap
  itk::ImageRegionConstIterator<TImage> nodeAIterator(this->Image, labelA_OverlapRegion);
  itk::ImageRegionConstIterator<TImage> nodeBIterator(this->Image, labelB_OverlapRegion);

  float ssd = 0.0f; // sum of squared differences
  while(!nodeAIterator.IsAtEnd())
    {
    // Get the value of the current pixel
    typename TImage::PixelType pixelA = nodeAIterator.Get();
    typename TImage::PixelType pixelB = nodeBIterator.Get();
    for(unsigned int i = 0; i < pixelA.GetSize(); ++i)
    {
      ssd += (pixelA[i] - pixelB[i]) * (pixelA[i] - pixelB[i]);
    }

    ++nodeAIterator;
    ++nodeBIterator;
    }

  return ssd/static_cast<float>(overlap.GetNumberOfPixels());
}

template <typename TImage>
float StructurePropagationDPConstrained<TImage>::UnaryEnergy(const LabelType& label, const NodeType& node)
{
//   assert(this->Image->GetLargestPossibleRegion().GetSize()[0] != 0);
//   assert(this->MaskImage->GetLargestPossibleRegion().GetSize()[0] != 0);
  assert(this->Image);
  assert(this->MaskImage);
  assert(this->MaskImage->IsValid(label));

  if(ITKHelpers::IndexDistance(ITKHelpers::GetRegionCenter(label), ITKHelpers::GetRegionCenter(node)) > this->MaxNodeLabelDistance)
  {
    return std::numeric_limits<float>::infinity();
  }
  // sum of the squared differences in the region of the patch which overlaps the source region

  if(this->MaskImage->IsHole(node))
  {
    return 0.0f;
  }

  itk::ImageRegionConstIterator<TImage> nodeIterator(this->Image, node);
  itk::ImageRegionConstIteratorWithIndex<Mask> maskIterator(this->MaskImage, node);
  itk::ImageRegionConstIterator<TImage> labelIterator(this->Image, label);

  float ssd = 0.0f; // sum of squared differences
  unsigned int numberOfPixelsUsed = 0;
  while(!nodeIterator.IsAtEnd())
    {
    // this is the region we want (assuming mask is black in the source region)
    if(this->MaskImage->IsValid(maskIterator.GetIndex()))
      {
      typename TImage::PixelType labelPixel = labelIterator.Get();
      typename TImage::PixelType nodePixel = nodeIterator.Get();
      for(unsigned int i = 0; i < labelPixel.GetSize(); ++i)
      {
        ssd += (labelPixel[i] - nodePixel[i]) * (labelPixel[i] - nodePixel[i]);
      }
      numberOfPixelsUsed++;
      }

    ++nodeIterator;
    ++maskIterator;
    ++labelIterator;
    }

  if(numberOfPixelsUsed == 0)
    {
    std::stringstream ss;
    ss << "Unary weight computation used 0 pixels! Node: " << node;
    throw std::runtime_error(ss.str());
    return 0.0f;
    }
  else
    {
    return this->UnaryBinaryWeight * ssd/static_cast<float>(numberOfPixelsUsed);
    }
}

#endif
