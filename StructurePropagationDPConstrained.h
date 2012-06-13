#ifndef StructurePropagationDPConstrained_h
#define StructurePropagationDPConstrained_h

// Submodules
#include "DynamicProgramming/DynamicProgramming.h"
#include "Mask/Mask.h"

// ITK
#include "itkImageRegion.h"

template <typename TImage>
class StructurePropagationDPConstrained : public DynamicProgramming<itk::ImageRegion<2>, itk::ImageRegion<2> >
{
public:

  StructurePropagationDPConstrained();

  typedef itk::ImageRegion<2> LabelType;
  typedef itk::ImageRegion<2> NodeType;

  void SetImage(TImage* const image);
  void SetMask(Mask* const mask);

  void SetMaxNodeLabelDistance(const unsigned int maxNodeLabelDistance);

  void SetUnaryBinaryWeight(const float unaryBinaryWeight);

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

  /** The maximum distance at which a label is allowed to be selected to fill a node. */
  unsigned int MaxNodeLabelDistance;

  /** The weight between unary and binary terms. */
  float UnaryBinaryWeight;
};

#include "StructurePropagationDPConstrained.hpp"

#endif
