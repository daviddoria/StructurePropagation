
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