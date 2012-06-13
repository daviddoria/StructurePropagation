#ifndef StructurePropagationDP_h
#define StructurePropagationDP_h

// Submodules
#include "DynamicProgramming/DynamicProgramming.h"
#include "Mask/Mask.h"

// ITK
#include "itkImageRegion.h"

template <typename TImage>
class StructurePropagationDP : public DynamicProgramming<itk::ImageRegion<2>, itk::ImageRegion<2> >
{
public:

  StructurePropagationDP();

  typedef itk::ImageRegion<2> LabelType;
  typedef itk::ImageRegion<2> NodeType;

  void SetImage(TImage* const image);
  void SetMask(Mask* const mask);

private:
  float BinaryEnergy(const LabelType& labelA, const NodeType& nodeA,
                     const LabelType& labelB, const NodeType& nodeB);

  float UnaryEnergy(const LabelType& label, const NodeType& node);

  /** The image from which to pull the patches. */
  //typename TImage::Pointer Image;
  TImage* Image;

  /** The mask from which to pull the patches. */
  //typename Mask::Pointer MaskImage;
  Mask* MaskImage;
};

#include "StructurePropagationDP.hpp"

#endif
