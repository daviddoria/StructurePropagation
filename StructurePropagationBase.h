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

#ifndef StructurePropagationBase_H
#define StructurePropagationBase_H

#include "Types.h"

class StructurePropagationBase
{
public:
  StructurePropagationBase();

  void ComputePatchRegions();
  void ComputeSourcePatchRegions();
  void ComputeTargetPatchRegions();

  void SetMask(UnsignedCharScalarImageType::Pointer mask);
  void SetPropagationLine(std::vector<itk::Index<2> > propagationLine);
  void SetPatchRadius(unsigned int radius);

  std::vector<itk::ImageRegion<2> > GetSourcePatchRegions();
  std::vector<itk::ImageRegion<2> > GetTargetPatchRegions();

  virtual void PropagateStructure() = 0;

protected:
  void ExtractRegionsOfPropagationPath();

  // Input data
  UnsignedCharScalarImageType::Pointer Mask;
  unsigned int PatchRadius;
  std::vector<itk::Index<2> > PropagationLine;

  // Derived data
  UnsignedCharScalarImageType::Pointer PropagationLineImage;

  std::vector<itk::ImageRegion<2> > SourcePatchRegions;
  std::vector<itk::ImageRegion<2> > TargetPatchRegions;

  std::vector<UnsignedCharScalarImageType::Pointer> TargetStrokePaths;
  std::vector<UnsignedCharScalarImageType::Pointer> SourceStrokePaths;

  // Keep track of the center of the nodes (why is this not itk::Index?)
  std::vector<itk::Point<int, 2> > NodeLocations;
};

#endif