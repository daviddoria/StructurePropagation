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

// ITK
#include "itkBinaryNotImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkContourMeanDistanceImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIterator.h"
#include "itkMaskImageFilter.h"
#include "itkPasteImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"

// Custom
#include "Helpers.h"

// OpenGM
#include <opengm/inference/treereweightedbeliefpropagation.hxx>
#include <opengm/inference/beliefpropagation.hxx>


template<class BP>
class BeliefPropagationProtocolVisitor {
public:
  typedef BP bp_type;
  typedef typename bp_type::value_type value_type;

  BeliefPropagationProtocolVisitor() : iteration(0){}

  void operator()(const bp_type& bp)
  {
    std::vector<size_t> state;
    bp.arg(state);

    value_type value = bp.graphicalModel().evaluate(state);

    value_type distance = bp.convergence();

    std::cout << "iteration " << iteration
        << ": energy=" <<  value
        << ", distance=" << distance
        << std::endl;
    iteration++;
  }

private:
  unsigned int iteration;

};

template <typename TImage>
StructurePropagation<TImage>::StructurePropagation()
{
  // Initializations
  this->Image = TImage::New();
  this->Mask = UnsignedCharScalarImageType::New();
  this->PropagationLineImage = UnsignedCharScalarImageType::New();
  this->PropagationLineSourceImage = UnsignedCharScalarImageType::New();
  this->PropagationLineTargetImage = UnsignedCharScalarImageType::New();

  // Default Values
  this->PatchRadius = 10;
  this->StructureWeight = 50.0;
  this->BoundaryMatchWeight = 2.0;
  this->SourceLineWidth = 3;
}

template <typename TImage>
typename TImage::Pointer StructurePropagation<TImage>::GetOutputImage()
{
  return this->OutputImage;
}

template <typename TImage>
void StructurePropagation<TImage>::SetImage(typename TImage::Pointer image)
{
  //this->Image = image;
  this->Image->Graft(image);
}

template <typename TImage>
void StructurePropagation<TImage>::SetSourceLineWidth(unsigned int width)
{
  this->SourceLineWidth = width;
}


template <typename TImage>
double StructurePropagation<TImage>::StructureCost(int node, int label)
{
  //std::cout << "There are " << this->TargetStrokePaths.size() << " TargetStrokePaths." << std::endl;

  //std::cout << "this->TargetStrokePaths[node] has " << Helpers::CountNonZeroPixels(this->TargetStrokePaths[node]) << " non zero pixels." << std::endl;
  //std::cout << "this->SourceStrokePaths[label] has " << Helpers::CountNonZeroPixels(this->SourceStrokePaths[label]) << " non zero pixels." << std::endl;

  assert(Helpers::CountNonZeroPixels(this->TargetStrokePaths[node]) > 0);
  assert(Helpers::CountNonZeroPixels(this->SourceStrokePaths[label]) > 0);

  typedef itk::ContourMeanDistanceImageFilter <UnsignedCharScalarImageType, UnsignedCharScalarImageType> ContourMeanDistanceImageFilterType;
  ContourMeanDistanceImageFilterType::Pointer contourMeanDistanceImageFilter =
    ContourMeanDistanceImageFilterType::New();
  contourMeanDistanceImageFilter->SetInput1(this->TargetStrokePaths[node]);
  contourMeanDistanceImageFilter->SetInput2(this->SourceStrokePaths[label]);
  contourMeanDistanceImageFilter->Update();

  return contourMeanDistanceImageFilter->GetMeanDistance();
}

template <typename TImage>
double StructurePropagation<TImage>::UnaryCost(int node, int label)
{
  double structureCost = StructureCost(node, label);
  double completionCost = CompletionCost(node, label);

  //std::cout << "Structure cost: " << structureCost << " " << "Completion cost: " << completionCost << std::endl;

  return this->StructureWeight * structureCost + this->BoundaryMatchWeight * completionCost;
}


template <typename TImage>
double StructurePropagation<TImage>::CompletionCost(int node, int label)
{
  // "sum of the normalized squared differences in the region of the patch which overlaps the source region"
  itk::ImageRegionConstIterator<TImage> imageIterator(this->Image, this->TargetPatchRegions[node]);
  itk::ImageRegionConstIterator<UnsignedCharScalarImageType> maskIterator(this->Mask,this->TargetPatchRegions[node]);

  itk::ImageRegionConstIterator<TImage> patchIterator(this->Image, this->SourcePatchRegions[label]);

  double ssd = 0; // sum of squared differences
  unsigned int numberOfPixels = 0;
  while(!imageIterator.IsAtEnd())
    {
    if(maskIterator.Get() == 0) // this is the region we want (assuming mask is black in the source region)
      {
      // Get the value of the current pixel
      typename TImage::PixelType pixel1 = imageIterator.Get();
      typename TImage::PixelType pixel2 = patchIterator.Get();
      ssd += Helpers::SquaredDifference(pixel1, pixel2);
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
    return ssd/static_cast<double>(numberOfPixels);
    }
}

template <typename TImage>
double StructurePropagation<TImage>::BinaryCost(int node1, int node2, int label1, int label2)
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
  itk::ImageRegionConstIterator<TImage> patch1Iterator(this->Image, overlapRegion1);
  itk::ImageRegionConstIterator<TImage> patch2Iterator(this->Image, overlapRegion2);

  double ssd = 0; // sum of squared differences
  while(!patch1Iterator.IsAtEnd())
    {
    // Get the value of the current pixel
    typename TImage::PixelType pixel1 = patch1Iterator.Get();
    typename TImage::PixelType pixel2 = patch2Iterator.Get();
    ssd += Helpers::SquaredDifference(pixel1, pixel2);

    ++patch1Iterator;
    ++patch2Iterator;
    }

  /*
  std::cout << "Binary cost for nodes " << node1 << " and " << node2
            << " with labels " << label1 << " and " << label2 << " is " << ssd << std::endl;
  */

  if(ssd <= 0)
    {
    return 0.0;
    }
  else
    {
    return ssd/static_cast<double>(overlapRegion1.GetNumberOfPixels());
    }
}


template <typename TImage>
void StructurePropagation<TImage>::CreateUnaryFactors(GraphicalModel &gm, Space &space)
{
  std::cout << "Creating unary factors..." << std::endl;

  unsigned int numberOfNodes = this->TargetPatchRegions.size();
  unsigned int numberOfLabels = this->SourcePatchRegions.size();

  //std::cout << "There are " << numberOfNodes << " nodes." << std::endl;
  //std::cout << "There are " << numberOfLabels << " labels." << std::endl;
  unsigned int numberOfUnaryFactors = 0;
  for(unsigned int node = 0; node < numberOfNodes; node++)
    {
    std::vector<size_t> unaryIndices(1, node);
    Factor unaryFactor(space, unaryIndices.begin(), unaryIndices.end());
    for(unsigned int label = 0; label < numberOfLabels; label++)
      {
      double cost = UnaryCost(node, label);
      unaryFactor(label) = cost;
      numberOfUnaryFactors++;
      //std::cout << "The cost of assigning node " << node << " = " << label << " is " << unaryFactor(label) << std::endl;
      }

    gm.addFactor(unaryFactor);
    }
  std::cout << "Finished creating unary factors! (There were " << numberOfUnaryFactors << ")" << std::endl;
}

template <typename TImage>
void StructurePropagation<TImage>::ClearEverything()
{
  this->SourcePatchRegions.clear();
  this->TargetPatchRegions.clear();
  this->NodeLocations.clear();
  this->Edges.clear();
}

template <typename TImage>
void StructurePropagation<TImage>::WriteTargetPatches()
{
  typename TImage::Pointer targetImage = TImage::New();
  targetImage->SetRegions(this->Image->GetLargestPossibleRegion());
  targetImage->Allocate();
  targetImage->FillBuffer(itk::NumericTraits<typename TImage::PixelType>::Zero);

  for(unsigned int i = 0; i < this->TargetPatchRegions.size(); i++)
    {
    Helpers::CopyPatchIntoImage<TImage>(this->Image, targetImage, this->TargetPatchRegions[i], this->TargetPatchRegions[i]);
    }
  Helpers::CastAndWriteImage<TImage>(targetImage, "TargetPatches.png");
}

template <typename TImage>
void StructurePropagation<TImage>::WriteSourcePatches()
{
  typename TImage::Pointer targetImage = TImage::New();
  targetImage->SetRegions(this->Image->GetLargestPossibleRegion());
  targetImage->Allocate();
  targetImage->FillBuffer(itk::NumericTraits<typename TImage::PixelType>::Zero);

  for(unsigned int i = 0; i < this->SourcePatchRegions.size(); i++)
    {
    Helpers::CopyPatchIntoImage<TImage>(this->Image, targetImage, this->SourcePatchRegions[i], this->SourcePatchRegions[i]);
    }
  Helpers::CastAndWriteImage<TImage>(targetImage, "SourcePatches.png");
}

template <typename TImage>
void StructurePropagation<TImage>::PropagateStructure()
{
  // Find intersections of target line with hole boundary and add them as additional intersections (so nodes are required to be placed there)
  for(std::set<itk::Index<2>, IndexComparison>::iterator iterator = this->PropagationLine.begin();
      iterator != this->PropagationLine.end(); iterator++)
    {
    if(this->Mask->GetPixel(*iterator)) // the pixel is in the target region
      {
      if(Helpers::HasNeighborWithValue<UnsignedCharScalarImageType>(this->Mask, *iterator, 0))
        {
        this->PropagationLineIntersections.insert(*iterator);
        std::cout << "Inserted boundary intersection." << std::endl;
        }
      }
    }

  // The order of the following functions is important
  CreateLineImages();
  CreateEdges();
  ComputePatchRegions();
  ExtractRegionsOfPropagationPath();

  WriteSourcePatches();
  WriteTargetPatches();

  Helpers::WritePixels(this->Image->GetLargestPossibleRegion().GetSize(), this->NodeLocations, "Nodes.png");

  // Create the graphical model
  unsigned int numberOfNodes = this->TargetPatchRegions.size();
  unsigned int numberOfLabels = this->SourcePatchRegions.size();
  std::cout << "PropagateStructure: There are " << numberOfNodes << " nodes." << std::endl;
  std::cout << "PropagateStructure: There are " << numberOfLabels << " labels." << std::endl;
  std::cout << "There are " << Helpers::CountNonZeroPixels(this->PropagationLineTargetImage) << " pixels in target lines.";
  std::cout << "The ratio of pixels to nodes is " << static_cast<double>(Helpers::CountNonZeroPixels(this->PropagationLineTargetImage))/
                                                      static_cast<double>(numberOfNodes);

  std::vector<size_t> nodes(numberOfNodes, numberOfLabels);
  Space space(nodes.begin(), nodes.end());
  GraphicalModel gm;

  CreateUnaryFactors(gm, space);
  CreateBinaryFactors(gm, space);

  ////////////// Optimize ///////////////
  typedef opengm::TreeReweightedBeliefPropagation<GraphicalModel, opengm::Minimizer, opengm::MaxDistance> BP;
  //typedef opengm::BeliefPropagation<GraphicalModel, opengm::Minimizer, opengm::MaxDistance> BP;
  std::cout << "Setting up tree-reweighted belief propagation... " << std::endl;
  BP::Parameter para;
  para.maximumNumberOfSteps_ = 100;

  BP bp(gm, para);
  BeliefPropagationProtocolVisitor<BP> visitor;
  bp.infer(visitor);

  std::vector<size_t> result;
  bp.arg(result);

  // Fill in the image with the best patches
  typedef itk::ImageDuplicator< TImage > ImageDuplicatorType;
  typename ImageDuplicatorType::Pointer duplicator = ImageDuplicatorType::New();
  duplicator->SetInputImage(this->Image);
  duplicator->Update();

  this->OutputImage = duplicator->GetOutput();

  for(unsigned int i = 0; i < numberOfNodes; i++)
    {
    Helpers::CopySelfPatchIntoTargetRegion<TImage>(this->OutputImage, this->Mask, this->SourcePatchRegions[result[i]],
                                           Helpers::GetCenteredRegionRadius(this->NodeLocations[i], this->PatchRadius));
    }

  //Helpers::CastAndWriteImage<TImage>(this->OutputImage, "result.png");

  std::cout << "Finished propagating structure!" << std::endl;
}



template <typename TImage>
void StructurePropagation<TImage>::CreateBinaryFactors(GraphicalModel &gm, Space &space)
{
  std::cout << "Creating binary factors..." << std::endl;

  unsigned int numberOfLabels = this->SourcePatchRegions.size();
  //std::cout << "There are " << numberOfLabels << " labels." << std::endl;
  //std::cout << "There are " << this->Edges.size() << " edges." << std::endl;
  //std::cout << "There are " << this->NodeLocations.size() << " nodes." << std::endl;

  unsigned int numberOfBinaryFactors = 0;
  for(unsigned int edgeId = 0; edgeId < this->Edges.size(); edgeId++)
    {
    std::pair<unsigned int, unsigned int> edge = this->Edges[edgeId];
    //std::cout << "Creating binary factor between " << edge.first << " and " << edge.second << std::endl;

    std::vector<size_t> binaryIndices(2);
    binaryIndices[0] = edge.first;
    binaryIndices[1] = edge.second;
    Factor binaryFactor(space, binaryIndices.begin(), binaryIndices.end());
    for(unsigned int labelA = 0; labelA < numberOfLabels; labelA++)
      {
      for(unsigned int labelB = 0; labelB < numberOfLabels; labelB++)
        {
        binaryFactor(labelA,labelB) = BinaryCost(edge.first, edge.second, labelA, labelB);
        numberOfBinaryFactors++;
        //std::cout << "The cost of assigning node " << node << " = " << labelA << " and node+1 " << node+1 << " = " << labelB << " is " << binaryFactor(labelA, labelB) << std::endl;
        }
      }
    gm.addFactor(binaryFactor);
    }

  std::cout << "Finished creating binary factors! (There are " << numberOfBinaryFactors << ")" << std::endl;
}

template <typename TImage>
void StructurePropagation<TImage>::SetMask(UnsignedCharScalarImageType::Pointer mask)
{
  this->Mask->Graft(mask);
}

template <typename TImage>
void StructurePropagation<TImage>::SetPropagationLine(std::set<itk::Index<2>, IndexComparison > propagationLine)
{
  this->PropagationLine = propagationLine;
}


template <typename TImage>
void StructurePropagation<TImage>::SetPropagationLineIntersections(std::set<itk::Index<2>, IndexComparison > intersections)
{
  this->PropagationLineIntersections = intersections;
}

template <typename TImage>
void StructurePropagation<TImage>::SetPatchRadius(unsigned int radius)
{
  this->PatchRadius = radius;
}

template <typename TImage>
void StructurePropagation<TImage>::CreateLineImages()
{
  Helpers::IndicesToBinaryImage(this->PropagationLine, this->Mask->GetLargestPossibleRegion(), this->PropagationLineImage);

  // Extract the part of the lines in the target region
  typedef itk::MaskImageFilter< UnsignedCharScalarImageType, UnsignedCharScalarImageType> MaskFilterType;
  MaskFilterType::Pointer targetMaskFilter = MaskFilterType::New();
  targetMaskFilter->SetInput(this->PropagationLineImage);
  targetMaskFilter->SetMaskImage(this->Mask);
  targetMaskFilter->Update();
  this->PropagationLineTargetImage->Graft(targetMaskFilter->GetOutput());

  Helpers::WriteImage<UnsignedCharScalarImageType>(this->PropagationLineTargetImage, "TargetLine.png");

  // Extract the part of the lines in the source region
  // Invert the mask
  typedef itk::BinaryNotImageFilter <UnsignedCharScalarImageType> BinaryNotImageFilterType;

  BinaryNotImageFilterType::Pointer binaryNotFilter = BinaryNotImageFilterType::New();
  binaryNotFilter->SetInput(this->Mask);
  binaryNotFilter->Update();

  MaskFilterType::Pointer sourceMaskFilter = MaskFilterType::New();
  sourceMaskFilter->SetInput(this->PropagationLineImage);
  sourceMaskFilter->SetMaskImage(binaryNotFilter->GetOutput());
  sourceMaskFilter->Update();
  this->PropagationLineSourceImage->Graft(sourceMaskFilter->GetOutput());

  Helpers::WriteImage<UnsignedCharScalarImageType>(this->PropagationLineSourceImage, "SourceLine.png");

}

template <typename TImage>
void StructurePropagation<TImage>::ComputeSourcePatchRegions()
{
  // Start fresh
  this->SourcePatchRegions.clear();

  UnsignedCharScalarImageType::Pointer dilatedLineImage = UnsignedCharScalarImageType::New();
  Helpers::GetDilatedImage(this->PropagationLineSourceImage, dilatedLineImage, this->SourceLineWidth/2);
  Helpers::WriteImage<UnsignedCharScalarImageType>(dilatedLineImage, "DilatedSourceLines.png");

  // The dilated source line image will bleed into the target region,
  // but this is ok because we ensure that only patches completely in the sourcce region are kept

  std::vector<itk::Index<2> > dilatedLinePixelList = Helpers::BinaryImageToPixelList(dilatedLineImage);

  itk::Size<2> patchSize;
  patchSize.Fill(this->PatchRadius*2 + 1);
  //std::cout << "ComputeSourcePatchRegions: patchRadius " << this->PatchRadius << std::endl;
  //std::cout << "ComputeSourcePatchRegions: patchSize " << patchSize << std::endl;

  for(unsigned int i = 0; i < dilatedLinePixelList.size(); i++)
    {
    // Create a region centered on the pixels
    itk::Index<2> corner = Helpers::GetCorner(dilatedLinePixelList[i], this->PatchRadius);
    itk::ImageRegion<2> region(corner, patchSize);

    // If the region is valid, add it to the list of source patches
    if(Helpers::IsValidPatch(region, this->Mask))
      {
      this->SourcePatchRegions.push_back(region);
      }
    }
}

template <typename TImage>
void StructurePropagation<TImage>::ComputeTargetPatchRegions()
{
  // Start fresh
  this->TargetPatchRegions.clear();

  itk::Size<2> patchSize;
  patchSize.Fill(this->PatchRadius*2 + 1);

  for(unsigned int i = 0; i < this->NodeLocations.size(); i++)
    {
    // Create a region centered on the pixels
    itk::Index<2> corner = Helpers::GetCorner(this->NodeLocations[i], this->PatchRadius);
    itk::ImageRegion<2> region(corner, patchSize);

    // Keep patches which intersect the target region
    if(Helpers::IsIntersectTargetRegion(region, this->Mask));
      {
      this->TargetPatchRegions.push_back(region);
      }
    }
}

template <typename TImage>
std::vector<itk::ImageRegion<2> > StructurePropagation<TImage>::GetSourcePatchRegions()
{
  return this->SourcePatchRegions;
}

template <typename TImage>
std::vector<itk::ImageRegion<2> > StructurePropagation<TImage>::GetTargetPatchRegions()
{
  return this->TargetPatchRegions;
}

template <typename TImage>
void StructurePropagation<TImage>::ExtractRegionsOfPropagationPath()
{
  // This function extracts the patches of the line so they can be used in the StructureCost computation

  // Output the line for testing
  //Helpers::WriteImage<UnsignedCharScalarImageType>(this->PropagationLineImage, "PropagationLine_ExtractRegionsOfPropagationPath.png");

  // Extract the regions of the propagation path
  this->SourceStrokePaths.clear();
  typedef itk::RegionOfInterestImageFilter< UnsignedCharScalarImageType, UnsignedCharScalarImageType > ROIFilterType;

  for(unsigned int i = 0; i < this->SourcePatchRegions.size(); i++)
    {
    ROIFilterType::Pointer roiFilter = ROIFilterType::New();
    roiFilter->SetRegionOfInterest(this->SourcePatchRegions[i]);
    roiFilter->SetInput(this->PropagationLineImage);
    roiFilter->Update();
    this->SourceStrokePaths.push_back(roiFilter->GetOutput());
    }

  this->TargetStrokePaths.clear();
  for(unsigned int i = 0; i < this->TargetPatchRegions.size(); i++)
    {
    ROIFilterType::Pointer roiFilter = ROIFilterType::New();
    roiFilter->SetRegionOfInterest(this->TargetPatchRegions[i]);
    roiFilter->SetInput(this->PropagationLineImage);
    roiFilter->Update();
    this->TargetStrokePaths.push_back(roiFilter->GetOutput());
    }

  std::cout << "There are " << this->TargetStrokePaths.size() << " TargetStrokePaths." << std::endl;
  std::cout << "There are " << this->SourceStrokePaths.size() << " SourceStrokePaths." << std::endl;
}

template <typename TImage>
void StructurePropagation<TImage>::ComputePatchRegions()
{
  this->ComputeSourcePatchRegions();
  this->ComputeTargetPatchRegions();
}


template <typename TImage>
void StructurePropagation<TImage>::CreateEdges()
{
  // Create a list of edges between node ids (indices of this->NodeLocations)

  // Number the pixels of the target line
  IntScalarImageType::Pointer numberedPixels = IntScalarImageType::New();
  Helpers::NumberPixels(this->PropagationLineTargetImage, numberedPixels);
  Helpers::WriteColorMappedImage<IntScalarImageType>(numberedPixels, "NumberedPixels.png");
  Helpers::WriteImage<IntScalarImageType>(numberedPixels, "NumberedPixels.mhd");

  IntScalarImageType::Pointer clusteredPixels = IntScalarImageType::New();
  Helpers::DeepCopy<IntScalarImageType>(numberedPixels, clusteredPixels);

  //std::vector<itk::Index<2> > intersections = Helpers::FindIntersections(this->PropagationLineTargetImage);
  std::cout << "There are " << this->PropagationLineIntersections.size() << " intersections." << std::endl;
  Helpers::GrowNumberedPixelClustersWithIntersections(clusteredPixels, clusteredPixels, this->PatchRadius, this->PropagationLineIntersections);

  //Helpers::GrowNumberedPixelClusters(clusteredPixels, clusteredPixels, this->PatchRadius);
  //Helpers::GrowNumberedPixelClusters(clusteredPixels, clusteredPixels, this->PatchRadius/2);

  /*
  std::cout << "CreateEdges: There are " << Helpers::CountNumberOfNonNegativeValues(clusteredPixels)
            << " nodes." << std::endl;
  */

  //this->NodeLocations = Helpers::GetLabelCenters(clusteredPixels);
  this->NodeLocations = Helpers::GetLabelCentersWithIntersections(clusteredPixels, this->PropagationLineIntersections);

  //std::cout << "CreateEdges: There are " << this->NodeLocations.size() << " nodes." << std::endl;

  // Create an iterator
  itk::Size<2> radius;
  radius.Fill(1);

  typedef itk::ConstNeighborhoodIterator<IntScalarImageType> NeighborhoodIteratorType;
  NeighborhoodIteratorType iterator(radius, clusteredPixels,
                                    clusteredPixels->GetLargestPossibleRegion());

  //std::vector<NeighborhoodIteratorType::OffsetType> neighbors;
  std::vector<NeighborhoodIteratorType::OffsetType> neighbors = Helpers::CreateNeighborOffsetsInRadius(radius);
  //std::cout << "There are " << neighbors.size() << " neighbors." << std::endl;

  while(!iterator.IsAtEnd())
    {
    // don't consider pixels which are not part of the line
    if(iterator.GetCenterPixel() < 0)
      {
      ++iterator;
      continue;
      }
    // Use the clustered pixels to determine the connectivity of the groups.
    // The edge's first and second labels will actually be used to map node cluster centers together (with the same connectivity as the original clusters)
    for(unsigned int i = 0; i < neighbors.size(); i++)
      {
      bool inBounds;
      iterator.GetPixel(neighbors[i], inBounds);
      if(inBounds)
        {
        if(iterator.GetPixel(neighbors[i]) >= 0)
          {
          std::pair<unsigned int, unsigned int> edge;
          edge.first = iterator.GetCenterPixel();
          edge.second = iterator.GetPixel(neighbors[i]);
          //std::cout << "Potential edge: " << edge.first << " " << edge.second << std::endl;

          if(edge.second > edge.first) // we only want one edge between a pair of nodes because we will use this to build a factor graph
            {
            //std::cout << "Creating edge between nodes " << edge.first << " and " << edge.second << std::endl;
            this->Edges.push_back(edge);
            }
          } // end if pixel is positive
        } // end if inbounds

      } // end neighbors loop

    ++iterator;
    } // end while

  //std::cout << "Created " << this->Edges.size() << " edges." << std::endl;

  // For testing only
  WriteEdges();
}

template <typename TImage>
void StructurePropagation<TImage>::WriteEdges()
{
  UnsignedCharScalarImageType::Pointer edgeImage = UnsignedCharScalarImageType::New();
  edgeImage->SetRegions(this->Mask->GetLargestPossibleRegion());
  edgeImage->Allocate();
  edgeImage->FillBuffer(0);

  //std::cout << "WriteEdges: There are " << this->Edges.size() << " edges to write." << std::endl;

  for(unsigned int i = 0; i < this->Edges.size(); i++)
    {
    Helpers::DrawLineInImage(edgeImage, this->NodeLocations[this->Edges[i].first], this->NodeLocations[this->Edges[i].second]);
    }

  Helpers::WriteImage<UnsignedCharScalarImageType>(edgeImage, "EdgeImage.png");
}
