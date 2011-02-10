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

/* Based on "Image Completion with Structure Propagation" by Jian Sun */

#ifndef StructurePropagation_H
#define StructurePropagation_H

#include "StructurePropagation.h"
#include "Types.h"

#include <opengm/explicitfactor.hxx>
#include <opengm/graphicalmodel.hxx>
#include <opengm/adder.hxx>
#include <opengm/inference/treereweightedbeliefpropagation.hxx>
#include <opengm/maxdistance.hxx>

typedef double Energy;
typedef opengm::DiscreteSpace Space;
typedef opengm::ExplicitFactor<Energy> Factor;
typedef opengm::GraphicalModel<Factor, opengm::Adder> GraphicalModel;

template <typename TImageType>
class StructurePropagation
{
public:
  StructurePropagation();

  void SetImage(typename TImageType::Pointer image);

  typename TImageType::Pointer GetOutputImage();

  void PropagateStructure();

  // Main functions
  void ComputePatchRegions(); // This function simply calls ComputeSourcePatchRegions and ComputeTargetPatchRegions
  void ComputeSourcePatchRegions();
  void ComputeTargetPatchRegions();

  // Mutators
  void SetMask(UnsignedCharScalarImageType::Pointer mask);
  void SetPropagationLine(std::vector<itk::Index<2> > propagationLine);
  void SetPatchRadius(unsigned int radius);

  // Accessors
  std::vector<itk::ImageRegion<2> > GetSourcePatchRegions();
  std::vector<itk::ImageRegion<2> > GetTargetPatchRegions();

private:

  double StructureCost(int node, int label);
  double UnaryCost(int node, int label);

  double CompletionCost(int node, int label);

  double BinaryCost(int node1, int node2, int label1, int label2);

  void CreateUnaryFactors(GraphicalModel &gm, Space &space);
  void CreateBinaryFactors(GraphicalModel &gm, Space &space);

  typename TImageType::Pointer Image;
  typename TImageType::Pointer OutputImage;

  // This function extracts the patches of the line so they can be used in the StructureCost computation
  void ExtractRegionsOfPropagationPath();
  std::vector<UnsignedCharScalarImageType::Pointer> TargetStrokePaths;
  std::vector<UnsignedCharScalarImageType::Pointer> SourceStrokePaths;

  // Input data
  UnsignedCharScalarImageType::Pointer Mask;
  unsigned int PatchRadius;
  std::vector<itk::Index<2> > PropagationLine;

  // Images of the lines drawn
  void CreateLineImages();
  UnsignedCharScalarImageType::Pointer PropagationLineImage;
  UnsignedCharScalarImageType::Pointer PropagationLineTargetImage;
  UnsignedCharScalarImageType::Pointer PropagationLineSourceImage;

  std::vector<itk::ImageRegion<2> > SourcePatchRegions;
  std::vector<itk::ImageRegion<2> > TargetPatchRegions;

  std::vector<itk::Index<2> > NodeLocations; // These indices correspond to those of TargetPatchRegions

  std::vector<std::pair<unsigned int, unsigned int> > Edges;
  void CreateEdges();
  void WriteEdges();
};

#include "StructurePropagation.txx"

#endif