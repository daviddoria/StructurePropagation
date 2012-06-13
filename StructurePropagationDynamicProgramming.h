#ifndef StructurePropagationDynamicProgramming_h
#define StructurePropagationDynamicProgramming_h

// Submodules
#include "DynamicProgramming/DynamicProgramming.h"

// ITK
#include "itkImageRegion.h"

class StructurePropagationDynamicProgramming : public DynamicProgramming<itk::ImageRegion<2> >
{
public:
  typedef itk::VectorImage<unsigned char, 2> ImageType;
  typedef itk::ImageRegion<2> LabelType;
  
  void SetImage(ImageType* const image);
  
private:
  float BinaryEnergy(const LabelType& labelA, const unsigned int nodeA,
                     const LabelType& labelB, const unsigned int nodeB);

  float UnaryEnergy(const LabelType& label, const unsigned int node);

  ImageType::Pointer Image;

  
};

#endif
