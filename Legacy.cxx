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