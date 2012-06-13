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

/* Based on "Image Completion with Structure Propagation" by Jian Sun SIGGRAPH 2005*/

#ifndef StructurePropagation_H
#define StructurePropagation_H

// Submodules
#include "DynamicProgramming/DynamicProgramming.h"
#include "Mask/Mask.h"

template <typename TImage>
class StructurePropagation
{
public:
  /** The type of the propagation line image. */
  typedef ITKHelpersTypes::UnsignedCharScalarImageType PropagationLineImageType;

  /** Constructor */
  StructurePropagation();

  /** Actually perform the structure propagation */
  void PropagateStructure();

  /** Clear everything. */
  void ClearEverything();

  /** This function simply calls ComputeSourcePatchRegions and ComputeTargetPatchRegions */
  void ComputePatchRegions();

  /** Compute the source/label regions from the propagation line image. */
  void ComputeSourcePatchRegions(const PropagationLineImageType* const propagationLineOutsideMask);

  /** Compute the target/node regions from the propagation line image. */
  void ComputeTargetPatchRegions(const PropagationLineImageType* const propagationLineInsideMask);

  /** Set the mask to indicate the hole region */
  void SetMask(Mask* mask);

  /** Set the radius of the patches to use. */
  void SetPatchRadius(const unsigned int radius);

  /** Set the image on which to operate. */
  void SetImage(TImage* const image);

  /** Get the regions considered as labels/source patches. */
  std::vector<itk::ImageRegion<2> > GetSourcePatchRegions();

  /** Get the regions considered as nodes/target patches. */
  std::vector<itk::ImageRegion<2> > GetTargetPatchRegions();

  /** Get the output of the algorithm (the inpainted image) */
  TImage* GetOutputImage();

  /** Set the image where non-zero values indicate the propagation line drawn by the user. */
  void SetPropagationLineImage(PropagationLineImageType* const propagationLineImage);

private:

  /** The input image */
  typename TImage::Pointer Image;

  /** The output image */
  typename TImage::Pointer OutputImage;

  /** The image where non-zero values indicate the propagation line drawn by the user. */
  PropagationLineImageType::Pointer PropagationLineImage;

  /** The mask indicating the hole in the image. */
  Mask::Pointer MaskImage;

  /** The radius of the patches to use. */
  unsigned int PatchRadius;

  /** The regions of the image corresponding to source/label patches. */
  std::vector<itk::ImageRegion<2> > SourcePatchRegions;

  /** The regions of the image corresponding to target/node patches. */
  std::vector<itk::ImageRegion<2> > TargetPatchRegions;

};

#include "StructurePropagation.hpp"

#endif