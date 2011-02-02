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

#include "itkCastImageFilter.h"
#include "itkContourMeanDistanceImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkImageFileWriter.h"
#include "itkPasteImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"

#include "Helpers.h"

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
StructurePropagation<TImageType>::StructurePropagation(typename TImageType::Pointer image)
{
  this->Image = TImageType::New();
  this->Image->Graft(image);
}

template <typename TImageType>
void StructurePropagation<TImageType>::SetImage(typename TImageType::Pointer image)
{
  //this->Image = image;
  this->Image->Graft(image);
}

template <typename TImageType>
double StructurePropagation<TImageType>::UnaryCost(int node, int label)
{
  double structureCost = StructureCost(node, label);
  double completionCost = CompletionCost(node, label);

  //std::cout << "Structure cost: " << structureCost << " " << "Completion cost: " << completionCost << std::endl;

  return structureCost + completionCost;
}

template <typename TImageType>
double StructurePropagation<TImageType>::StructureCost(int node, int label)
{
  //std::cout << "There are " << this->TargetStrokePaths.size() << " TargetStrokePaths." << std::endl;

  typedef itk::ContourMeanDistanceImageFilter <UnsignedCharScalarImageType, UnsignedCharScalarImageType> ContourMeanDistanceImageFilterType;

  //std::cout << "this->TargetStrokePaths[node] has " << Helpers::CountNonZeroPixels(this->TargetStrokePaths[node]) << " non zero pixels." << std::endl;
  //std::cout << "this->SourceStrokePaths[label] has " << Helpers::CountNonZeroPixels(this->SourceStrokePaths[label]) << " non zero pixels." << std::endl;

  assert(Helpers::CountNonZeroPixels(this->TargetStrokePaths[node]) > 0);
  assert(Helpers::CountNonZeroPixels(this->SourceStrokePaths[label]) > 0);

  typename ContourMeanDistanceImageFilterType::Pointer contourMeanDistanceImageFilter =
    ContourMeanDistanceImageFilterType::New();
  contourMeanDistanceImageFilter->SetInput1(this->TargetStrokePaths[node]);
  contourMeanDistanceImageFilter->SetInput2(this->SourceStrokePaths[label]);
  contourMeanDistanceImageFilter->Update();

  return contourMeanDistanceImageFilter->GetMeanDistance();
}

template <typename TImageType>
double StructurePropagation<TImageType>::CompletionCost(int node, int label)
{
  // "sum of the normalized squared differences in the region of the patch which overlaps the source region"
  itk::ImageRegionConstIterator<TImageType> imageIterator(this->Image, this->TargetPatchRegions[node]);
  itk::ImageRegionConstIterator<UnsignedCharScalarImageType> maskIterator(this->Mask,this->TargetPatchRegions[node]);

  itk::ImageRegionConstIterator<TImageType> patchIterator(this->Image, this->SourcePatchRegions[label]);

  double ssd = 0; // sum of squared differences
  double numberOfPixels = 0;
  while(!imageIterator.IsAtEnd())
    {
    if(maskIterator.Get() == 0) // this is the region we want (assuming mask is black in the source region)
      {
      // Get the value of the current pixel
      typename TImageType::PixelType pixel1 = imageIterator.Get();
      typename TImageType::PixelType pixel2 = patchIterator.Get();
      ssd += Helpers::difference(pixel1, pixel2);
      numberOfPixels++;
      }

    ++imageIterator;
    ++maskIterator;
    ++patchIterator;
    }
  if(numberOfPixels <= 0)
    {
    return 0;
    }
  else
    {
    return ssd/numberOfPixels;
    }
}

template <typename TImageType>
double StructurePropagation<TImageType>::BinaryCost(int node1, int node2, int label1, int label2)
{
  itk::ImageRegion<2> targetRegion1 = this->TargetPatchRegions[node1];
  itk::ImageRegion<2> targetRegion2 = this->TargetPatchRegions[node2];

  itk::ImageRegion<2> sourceRegion1 = this->SourcePatchRegions[label1];
  itk::ImageRegion<2> sourceRegion2 = this->SourcePatchRegions[label2];

  // Find the overlapping region
  itk::ImageRegion<2> overlap = targetRegion1;
  overlap.Crop(targetRegion2);

  // Find the offset of the corner of the overlap region from each of the patches
  itk::Offset<2> offset1 = overlap.GetIndex() - targetRegion1.GetIndex();
  itk::Offset<2> offset2 = overlap.GetIndex() - targetRegion2.GetIndex();

  itk::ImageRegion<2> overlapRegion1(sourceRegion1.GetIndex() + offset1, overlap.GetSize());
  itk::ImageRegion<2> overlapRegion2(sourceRegion2.GetIndex() + offset2, overlap.GetSize());

  // sum of the normalized squared differences in the region where the two patches overlap
  itk::ImageRegionConstIterator<TImageType> patch1Iterator(this->Image, overlapRegion1);
  itk::ImageRegionConstIterator<TImageType> patch2Iterator(this->Image, overlapRegion2);

  double ssd = 0; // sum of squared differences
  while(!patch1Iterator.IsAtEnd())
    {
    // Get the value of the current pixel
    typename TImageType::PixelType pixel1 = patch1Iterator.Get();
    typename TImageType::PixelType pixel2 = patch2Iterator.Get();
    //ssd += difference<typename TImageType::PixelType>(pixel1, pixel2);
    ssd += Helpers::difference(pixel1, pixel2);

    ++patch1Iterator;
    ++patch2Iterator;
    }

  std::cout << "Binary cost for nodes " << node1 << " and " << node2 << " with labels " << label1 << " and " << label2 << " is " << ssd << std::endl;

  if(ssd <= 0)
    {
    return 0.0;
    }
  else
    {
    return ssd/static_cast<double>(overlapRegion1.GetNumberOfPixels());
    }
}

template <typename TImageType>
void StructurePropagation<TImageType>::PropagateStructure()
{
  //std::cout << "TImageType has " << TImageType::PixelType::GetNumberOfComponents() << " components." << std::endl;

  this->PropagationLineImage->SetRegions(this->Mask->GetLargestPossibleRegion());
  this->PropagationLineImage->Allocate();
  Helpers::IndicesToBinaryImage(this->PropagationLine, this->PropagationLineImage);

  ComputePatchRegions();
  ExtractRegionsOfPropagationPath();

  // Create the graphical model
  unsigned int numberOfNodes = this->TargetPatchRegions.size();
  unsigned int numberOfLabels = this->SourcePatchRegions.size();

  std::vector<size_t> nodes(numberOfNodes, numberOfLabels);
  Space space(nodes.begin(), nodes.end());
  GraphicalModel gm;

  // Create all unary factors
  std::cout << "Creating unary factors..." << std::endl;
  //std::cout << "There are " << numberOfNodes << " nodes." << std::endl;
  //std::cout << "There are " << numberOfLabels << " labels." << std::endl;
  for(unsigned int node = 0; node < numberOfNodes; node++)
    {
    std::vector<size_t> unaryIndices(1, node);
    Factor unaryFactor(space, unaryIndices.begin(), unaryIndices.end());
    for(unsigned int label = 0; label < numberOfLabels; label++)
      {
      double cost = UnaryCost(node, label);
      unaryFactor(label) = cost;
      //std::cout << "The cost of assigning node " << node << " = " << label << " is " << unaryFactor(label) << std::endl;
      }

    gm.addFactor(unaryFactor);
    }

  // Create all binary factors
  std::cout << "Creating binary factors..." << std::endl;
  for(unsigned int node = 0; node < numberOfNodes - 1; node++) // we loop until numberOfNodes-1 because we use the current node and the current node + 1 as the two nodes in the factor
    {
    std::vector<size_t> binaryIndices(2);
    binaryIndices[0] = node;
    binaryIndices[1] = node+1;
    Factor binaryFactor(space, binaryIndices.begin(), binaryIndices.end());
    for(unsigned int labelA = 0; labelA < numberOfLabels; labelA++)
      {
      for(unsigned int labelB = 0; labelB < numberOfLabels; labelB++)
        {
        binaryFactor(labelA,labelB) = BinaryCost(node, node+1, labelA, labelB);
        //std::cout << "The cost of assigning node " << node << " = " << labelA << " and node+1 " << node+1 << " = " << labelB << " is " << binaryFactor(labelA, labelB) << std::endl;
        }
      }
    gm.addFactor(binaryFactor);
    }

  ////////////// Optimize ///////////////
  typedef opengm::TreeReweightedBeliefPropagation<GraphicalModel, opengm::Minimizer, opengm::MaxDistance> TRBP;
  std::cout << "Setting up tree-reweighted belief propagation... " << std::endl;
  TRBP::Parameter para;
  para.maximumNumberOfSteps_ = 100;

  TRBP trbp(gm, para);
  trbp.infer();

  std::vector<size_t> result;
  trbp.arg(result);

  std::cout << "Best labeling:" << std::endl;
  for(unsigned int i = 0; i < result.size(); i++)
    {
    std::cout << "Node " << i << " label " << result[i] << std::endl;
    }

  // Fill in the image with the best patches
  typedef itk::ImageDuplicator< TImageType > ImageDuplicatorType;
  typename ImageDuplicatorType::Pointer duplicator = ImageDuplicatorType::New();
  duplicator->SetInputImage(this->Image);
  duplicator->Update();

  typename TImageType::Pointer outputImage = duplicator->GetOutput();

  typedef itk::PasteImageFilter <TImageType, TImageType >
    PasteImageFilterType;

  for(unsigned int i = 0; i < numberOfNodes; i++)
    {
    typename PasteImageFilterType::Pointer pasteFilter
      = PasteImageFilterType::New();
    pasteFilter->SetSourceImage(outputImage);
    pasteFilter->SetDestinationImage(outputImage);
    pasteFilter->SetSourceRegion(this->SourcePatchRegions[result[i]]);
    itk::Index<2> destinationIndex;
    destinationIndex[0] = this->NodeLocations[i][0] - this->PatchRadius;
    destinationIndex[1] = this->NodeLocations[i][1] - this->PatchRadius;
    pasteFilter->SetDestinationIndex(destinationIndex);
    pasteFilter->Update();

    outputImage->Graft(pasteFilter->GetOutput());
    }

  typedef itk::CovariantVector<unsigned char,TImageType::PixelType::Dimension> UnsignedCharVectorPixelType;
  typedef itk::Image<UnsignedCharVectorPixelType, 2> UnsignedCharVectorImageType;

  typedef itk::CastImageFilter< TImageType, UnsignedCharVectorImageType > CastFilterType;
  typename CastFilterType::Pointer castFilter = CastFilterType::New();
  castFilter->SetInput(outputImage);

  typedef  itk::ImageFileWriter< UnsignedCharVectorImageType > WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName("result.png");
  writer->SetInput(castFilter->GetOutput());
  writer->Update();

}
