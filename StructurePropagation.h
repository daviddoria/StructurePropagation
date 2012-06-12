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

template <typename TImageType>
class StructurePropagation
{
public:
  StructurePropagation();

  // Functions
  void PropagateStructure();
  void ClearEverything();

  // Main functions
  void ComputePatchRegions(); // This function simply calls ComputeSourcePatchRegions and ComputeTargetPatchRegions
  void ComputeSourcePatchRegions();
  void ComputeTargetPatchRegions();

  // Mutators
  void SetMask(Mask* mask);
  //void SetPropagationLine(std::vector<itk::Index<2> > propagationLine);
  void SetPropagationLine(std::set<itk::Index<2>, IndexComparison > propagationLine);
  void SetPropagationLineIntersections(std::set<itk::Index<2>, IndexComparison > intersections);
  void SetPatchRadius(unsigned int radius);
  void SetSourceLineWidth(unsigned int width);
  void SetImage(typename TImageType::Pointer image);

  // Accessors
  std::vector<itk::ImageRegion<2> > GetSourcePatchRegions();
  std::vector<itk::ImageRegion<2> > GetTargetPatchRegions();
  typename TImageType* GetOutputImage();

private:

  // Cost functions
  double UnaryCost(const int node, const int label);

  double CompletionCost(const int node, const int label);
  double BinaryCost(const int node1, const int node2, const int label1, const int label2);

  // The input and output images
  typename TImageType::Pointer Image;
  typename TImageType::Pointer OutputImage;

  // This function extracts the patches of the line so they can be used in the StructureCost computation
  void ExtractRegionsOfPropagationPath();
  std::vector<UnsignedCharScalarImageType::Pointer> TargetStrokePaths;
  std::vector<UnsignedCharScalarImageType::Pointer> SourceStrokePaths;

  // Data
  UnsignedCharScalarImageType::Pointer Mask;
  unsigned int PatchRadius;
  std::set<itk::Index<2>, IndexComparison > PropagationLine;
  std::set<itk::Index<2>, IndexComparison > PropagationLineIntersections;

  // Images of the propagation path.
  void CreateLineImages();
  UnsignedCharScalarImageType::Pointer PropagationLineImage;
  UnsignedCharScalarImageType::Pointer PropagationLineTargetImage;
  UnsignedCharScalarImageType::Pointer PropagationLineSourceImage;

  // The regions of the image corresponding to source and target patches.
  std::vector<itk::ImageRegion<2> > SourcePatchRegions;
  std::vector<itk::ImageRegion<2> > TargetPatchRegions;

  // These indices are the center of the patch associated with each node.
  // These indices correspond to those of TargetPatchRegions.
  std::vector<itk::Index<2> > NodeLocations;

  void WriteTargetPatches();
  void WriteSourcePatches();
};

#include "StructurePropagation.hpp"

#endif