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

void StructurePropagationBase::ComputeTargetPatchRegions()
{
  this->TargetPatchRegions.clear();

  // This function assumes that the stroke starts outside the masked region, enters the masked region, then exits the masked region
  bool enteredMaskedRegion = false;

  itk::Size<2> patchSize;
  patchSize.Fill(this->PatchRadius*2 + 1);

  for(unsigned int i = 0; i < this->PropagationLine.size(); i++)
    {
    // Create a region centered on the pixels
    itk::Index<2> corner = Helpers::GetCorner(this->PropagationLine[i], this->PatchRadius);
    itk::ImageRegion<2> region(corner, patchSize);

    // Keep patches which intersect the mask by finding the first one, then moving PatchRadius and continuing
    itk::ImageRegionConstIterator<UnsignedCharScalarImageType> imageIterator(this->Mask, region);

    bool inMaskedRegion = false;
    while(!imageIterator.IsAtEnd())
      {
      if(imageIterator.Get() > 0) // the mask is non-zero
        {
        inMaskedRegion = true;
        enteredMaskedRegion = true; // this is a flag that gets set once during the whole loop
        break; // stop checking the current patch
        }
      ++imageIterator;
      }

    // Add the region to the list of target patches
    if(inMaskedRegion)
      {
      this->TargetPatchRegions.push_back(region);
      i += this->PatchRadius; // skip the next few pixels
      /*
      itk::Point<int, 2> point;
      point[0] = this->PropagationLine[i][0];
      point[1] = this->PropagationLine[i][1];
      NodeLocations.push_back(point);
      */
      NodeLocations.push_back(this->PropagationLine[i]);

      }

    // Stop when we find the first patch that is not in the target region after we have gone through the hole
    if(!inMaskedRegion && enteredMaskedRegion)
      {
      break; // we are done
      }

    }
}


void ClusterNumberedPixels(IntScalarImageType::Pointer numberedPixels, IntScalarImageType::Pointer clusteredPixels)
{
  std::cout << "Starting ClusterNumberedPixels..." << std::endl;

  // Input: numberedPixels - an image where each valid pixel is an id >= 0
  // Output: clusteredPixels - a modified version of numberedPixels where the regions of each id has grown
  // Details: The function tries to expand the region associated with each label

  // Initialize by making the output the same as the input
  DeepCopy<IntScalarImageType>(numberedPixels, clusteredPixels);

  // This bool image will track which pixels have already been expanded or relabeled
  itk::Image<bool, 2>::Pointer visitedImage = itk::Image<bool, 2>::New();
  visitedImage->SetRegions(numberedPixels->GetLargestPossibleRegion());
  visitedImage->Allocate();
  visitedImage->FillBuffer(0);

  std::vector<itk::Index<2> > validPixels = GetNonNegativePixels(clusteredPixels);
  //std::cout << "There are " << validPixels.size() << " validPixels." << std::endl;

  for(unsigned int currentPixelId = 0; currentPixelId < validPixels.size(); currentPixelId++) // loop over line pixels
    {
    // Skip pixels which have already been visited
    if(visitedImage->GetPixel(validPixels[currentPixelId]))
      {
      continue;
      }
    visitedImage->SetPixel(validPixels[currentPixelId], 1); // mark as visited
    int currentValue = numberedPixels->GetPixel(validPixels[currentPixelId]);
    //std::cout << "currentValue: " << currentValue << std::endl;

    std::vector<itk::Index<2> > validNeighbors = GetNonNegativeValidNeighbors(clusteredPixels, validPixels[currentPixelId]);
    //std::cout << "There are " << validNeighbors.size() << " validNeighbors." << std::endl;

    for(unsigned int neighborPixelId = 0; neighborPixelId < validNeighbors.size(); neighborPixelId++) // loop over neighbors of the current pixel
      {
      // Skip pixels which have already been visited
      if(visitedImage->GetPixel(validNeighbors[neighborPixelId]))
        {
        continue;
        }

      // Don't mark neighbors with the same label as visited
      if(numberedPixels->GetPixel(validNeighbors[neighborPixelId]) == currentValue)
        {
        continue;
        }

      //std::cout << "Changing " << clusteredPixels->GetPixel(iterator.GetIndex() + neighbors[i])
        //        << " to " << iterator.GetCenterPixel() << std::endl;
      clusteredPixels->SetPixel(validNeighbors[neighborPixelId], currentValue); // expand label
      visitedImage->SetPixel(validNeighbors[neighborPixelId], 1); // mark as visited
      //numberOfChangedPixels++;
      } // end neighbor loop
    } // end line pixel loop

  //std::cout << "There were " << numberOfChangedPixels << " changed pixels." << std::endl;
  Helpers::WriteColorMappedImage<IntScalarImageType>(clusteredPixels, "ClusteredPixelsBeforeRelabel.png");

  RelabelSequential(clusteredPixels);

  std::cout << "Finished ClusterNumberedPixels!" << std::endl;
}


void ClusterHalfNumberedPixels(IntScalarImageType::Pointer numberedPixels, IntScalarImageType::Pointer clusteredPixels)
{
  std::cout << "Starting ClusterNumberedPixels..." << std::endl;

  // Input: numberedPixels - an image where each valid pixel is an id >= 0
  // Output: clusteredPixels - a modified version of numberedPixels where the regions of each id has grown
  // Details: The function expand the region associated with each label


  // Initialize by making the output the same as the input
  DeepCopy<IntScalarImageType>(numberedPixels, clusteredPixels);

  std::set<int> valuesSet = GetNonNegativeUniqueValues<IntScalarImageType>(clusteredPixels);
  std::vector<unsigned int> valuesVector(valuesSet.begin(), valuesSet.end());

  for(unsigned int i = 0; i < valuesVector.size(); i+=2)
    {
    ChangeValue<IntScalarImageType>(clusteredPixels,valuesVector[i], -2); // mark pixels to be filled with -2
    }

  std::vector<itk::Index<2> > pixelsToFill = FindPixelsWithValue<IntScalarImageType>(clusteredPixels, -2);
  std::vector<itk::Index<2> > validPixels = GetNonNegativePixels(clusteredPixels);

  for(unsigned int currentPixelId = 0; currentPixelId < pixelsToFill.size(); currentPixelId++) // loop over line pixels
    {
    unsigned int closestIndexId = FindClosestIndex(validPixels, pixelsToFill[currentPixelId]);
    clusteredPixels->SetPixel(pixelsToFill[currentPixelId], clusteredPixels->GetPixel(validPixels[closestIndexId]));
    }

  RelabelSequential(clusteredPixels);

  std::cout << "Finished ClusterNumberedPixels!" << std::endl;
}