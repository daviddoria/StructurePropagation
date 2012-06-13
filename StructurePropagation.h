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

  /** This function divides the propagation line into the part of the line that is inside the hole
    * and the part of the line that is outside the hole. It then calls ComputeSourceRegions and ComputeTargetRegions
    * on the corrsponding parts of the line. */
  void ComputeRegionsPropagationThroughHole();

  /** This function computes the source region in a dilated ring around the hole, and the target regions along the hole boundary. */
  void ComputeRegionsPropagationAroundHole();

  /** Set the regions considered as labels/source patches. */
  void SetSourceRegions(const std::vector<itk::ImageRegion<2> >& sourceRegions);

  /** Set the regions considered as nodes/target patches. */
  void SetTargetRegions(const std::vector<itk::ImageRegion<2> >& targetRegions);

  /** Compute the source/label regions from the propagation line image. */
  void ComputeSourceRegions(const PropagationLineImageType* const propagationLineOutsideMask);

  /** Compute the target/node regions from the propagation line image. */
  void ComputeTargetRegionsOpenContour(const PropagationLineImageType* const propagationLineInsideMask);

  /** Compute the target/node regions from the propagation line image. */
  void ComputeTargetRegionsClosedContour(const PropagationLineImageType* const propagationLineInsideMask);

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
  std::vector<itk::ImageRegion<2> > SourceRegions;

  /** The regions of the image corresponding to target/node patches. */
  std::vector<itk::ImageRegion<2> > TargetRegions;

};

#include "StructurePropagation.hpp"

#endif