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

  virtual void ExtractSourcePatches(std::vector<itk::Index<2> > pixels) = 0;
  virtual void ComputeTargetPatches(std::vector<itk::Index<2> > pixels) = 0;

  void SetMask(UnsignedCharScalarImageType::Pointer mask);

  void SetPatchRadius(unsigned int radius);

protected:
  UnsignedCharScalarImageType::Pointer Mask;
  unsigned int PatchRadius;

  std::vector<itk::ImageRegion<2> > PatchRegions;
};

#endif