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

#ifndef TYPES_H
#define TYPES_H

#include "itkImage.h"
#include "itkCovariantVector.h"

// All images are stored internally as float pixels
typedef itk::CovariantVector<float,3> ColorPixelType;
typedef itk::CovariantVector<float,1> GrayscalePixelType;

typedef itk::Image<ColorPixelType, 2> ColorImageType;
typedef itk::Image<GrayscalePixelType, 2> GrayscaleImageType;

//typedef itk::Image<GrayscalePixelType, 2> MaskImageType;

// For writing images, we need to first convert to actual unsigned char images
typedef itk::Image<unsigned char, 2> UnsignedCharScalarImageType;
typedef UnsignedCharScalarImageType MaskImageType;

typedef itk::Image<int, 2> IntScalarImageType;


#endif