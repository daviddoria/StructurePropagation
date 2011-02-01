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

  void ComputePatchRegions(std::vector<itk::Index<2> > pixels);
  void ComputeSourcePatchRegions(std::vector<itk::Index<2> > pixels);
  void ComputeTargetPatchRegions(std::vector<itk::Index<2> > pixels);

  void SetMask(UnsignedCharScalarImageType::Pointer mask);
  void GetMask(UnsignedCharScalarImageType::Pointer mask);

  void SetPropagationLine(UnsignedCharScalarImageType::Pointer propagationLine);

  void SetPatchRadius(unsigned int radius);

  std::vector<itk::ImageRegion<2> > GetSourcePatchRegions();
  std::vector<itk::ImageRegion<2> > GetTargetPatchRegions();

  virtual void PropagateStructure() = 0;

protected:
  void ExtractRegionsOfPropagationPath();

  UnsignedCharScalarImageType::Pointer Mask;
  UnsignedCharScalarImageType::Pointer PropagationLine;
  unsigned int PatchRadius;

  std::vector<itk::ImageRegion<2> > SourcePatchRegions;
  std::vector<itk::ImageRegion<2> > TargetPatchRegions;

  std::vector<UnsignedCharScalarImageType::Pointer> TargetStrokePaths;
  std::vector<UnsignedCharScalarImageType::Pointer> SourceStrokePaths;

  //std::vector<itk::Index<2> > NodeLocations;
  std::vector<itk::Point<int, 2> > NodeLocations;
};

#endif