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

void UpdateImage()
{

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
      void CopySelfPatchIntoTargetRegion(typename TImageType::Pointer image, MaskImageType::Pointer mask, itk::ImageRegion<2> sourceRegion, itk::ImageRegion<2> targetRegion)
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
}

void CreateBinaryFactors()
{

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
}


template <typename TImage>
void StructurePropagation<TImage>::CreateEdges()
{
  // Create a list of edges between node ids (indices of this->NodeLocations)

  // Number the pixels of the target line
  IntScalarImageType::Pointer numberedPixels = IntScalarImageType::New();
  numberedPixels->SetRegions(this->PropagationLineTargetImage->GetLargestPossibleRegion());
  numberedPixels->Allocate();
  numberedPixels->FillBuffer(-1);

  Helpers::NumberPixels(this->PropagationLineTargetImage, numberedPixels);

  //Helpers::WriteImage<IntScalarImageType>(numberedPixels, "NumberedPixels.mhd");
  Helpers::WriteColorMappedImage<IntScalarImageType>(numberedPixels, "NumberedPixels.png");

  // Create an iterator
  itk::Size<2> radius;
  radius.Fill(1);

  typedef itk::ConstNeighborhoodIterator<IntScalarImageType> NeighborhoodIteratorType;
  NeighborhoodIteratorType iterator(radius, numberedPixels,
                                    numberedPixels->GetLargestPossibleRegion());

  std::vector<NeighborhoodIteratorType::OffsetType> neighbors;

  for(int i = -1 * static_cast<int>(radius[0]); i <= static_cast<int>(radius[0]); i++)
    {
    for(int j = -1 * static_cast<int>(radius[1]); j <= static_cast<int>(radius[1]); j++)
      {
      if(!(i==0 && j==0))
        {
        NeighborhoodIteratorType::OffsetType offset = {{i,j}};
        neighbors.push_back(offset);
        } // end if
      } // end j loop
    } // end i loop

  //std::cout << "There are " << neighbors.size() << " neighbors." << std::endl;



  while(!iterator.IsAtEnd())
    {
    // don't consider pixels which are not part of the line
    if(iterator.GetCenterPixel() < 0)
      {
      ++iterator;
      continue;
      }

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

          if(edge.second > edge.first) // we only want one edge between a pair of nodes, because we will use this to build a factor graph
            {
            this->Edges.push_back(edge);
            }
          } // end if pixel is positive
        } // end if inbounds

      } // end neighbors loop

    ++iterator;
    } // end while

  std::cout << "Created " << this->Edges.size() << " edges." << std::endl;

  // For testing only
  WriteEdges();
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
  for(unsigned int i = 0; i < 2; i++)
    {
    Helpers::ClusterNumberedPixels(clusteredPixels, clusteredPixels);
    //Helpers::WriteColorMappedImage<IntScalarImageType>(clusteredPixels, "ClusteredPixels1.png");
    std::cout << "CreateEdges: There are " << Helpers::CountNumberOfNonNegativeValues(clusteredPixels)
              << " nodes after iteration " << i << "." << std::endl;
    }

  //Helpers::WriteImage<IntScalarImageType>(clusteredPixels, "ClusteredPixels.mhd");

  //// for testing only ////
  /*
  for(unsigned int i = 0; i < Helpers::CountNumberOfNonNegativeValues(clusteredPixels); i++)
    {
    unsigned int number = Helpers::CountPixelsWithValue<IntScalarImageType>(clusteredPixels, i);
    std::cout << "CreateEdges Test: Label " << i << " has " << number << " pixels." << std::endl;
    }
  */
  /////

  this->NodeLocations = Helpers::GetLabelCenters(clusteredPixels);

  //std::cout << "CreateEdges: There are " << this->NodeLocations.size() << " nodes." << std::endl;

  // Create an iterator
  itk::Size<2> radius;
  radius.Fill(1);

  typedef itk::ConstNeighborhoodIterator<IntScalarImageType> NeighborhoodIteratorType;
  NeighborhoodIteratorType iterator(radius, clusteredPixels,
                                    clusteredPixels->GetLargestPossibleRegion());

  //std::vector<NeighborhoodIteratorType::OffsetType> neighbors;
  std::vector<NeighborhoodIteratorType::OffsetType> neighbors = Helpers::CreateOffsetsInRadius(radius);
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