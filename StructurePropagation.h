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

#ifndef StructurePropagation_H
#define StructurePropagation_H

#include "StructurePropagationBase.h"
#include "Types.h"

template <typename TImageType>
class StructurePropagation : public StructurePropagationBase
{
public:
  StructurePropagation(typename TImageType::Pointer image);

  void SetImage(typename TImageType::Pointer image);

  void PropagateStructure();

private:
  typename TImageType::Pointer Image;

  typename TImageType::Pointer OutputImage;

  double UnaryCost(int node, int label);
  double StructureCost(int node, int label);
  double CompletionCost(int node, int label);

  double BinaryCost(int node1, int node2, int label1, int label2);
};

#include "StructurePropagation.txx"

#endif