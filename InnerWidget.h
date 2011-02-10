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

/* This is the main GUI class of this project. It is a QMainWindow
 * so that we can use a File menu. It contains an instance of our main functional
 * class ImageGraphCutBase and our custom scribble interactor style vtkGraphCutInteractorStyle.
 * It also contains a CProgressThread so that we can display a progress bar in marquee
 * mode during long computations.
*/

#ifndef INNERWIDGET_H
#define INNERWIDGET_H

// ITK
#include "itkCovariantVector.h"
#include "itkImage.h"
#include "itkImageRegionConstIteratorWithIndex.h"

// Custom
#include "StructurePropagation.h"
#include "vtkScribbleInteractorStyle.h"
//#include "ProgressThread.h"

// VTK
#include <vtkSmartPointer.h>

// Boost.Signals
#include <boost/signals2/signal.hpp>
#include <boost/bind.hpp>

#include <string>

// Forward declarations
class vtkImageActor;
class vtkRenderer;

#include "InnerWidgetObject.h"

template <typename TImage>
class InnerWidget : public InnerWidgetObject
{
public:
  InnerWidget(QWidget *parent = 0);

  // Use a QFileDialog to get a filename, then open the specified file as a greyscale or color image, depending on which type the user has specified through the file menu.
  void OpenFile();

  // Menu items
  //void actionFlip_Image_triggered() = 0;
  //void actionSave_Result_triggered() = 0;

  // Buttons
  void btnClearStrokes_clicked();
  void btnSaveStrokes_clicked();
  void btnPropagate_clicked();
  void btnLoadMask_clicked();

  // Radio buttons
  //void radDrawPropagationLine_clicked();
  //void radDrawHole_clicked();

  // These slots handle running the progress bar while the computations are done in a separate thread.
  //void StartProgressSlot();
  //void StopProgressSlot();

protected:

  void DisplayImage(typename TImage::Pointer image);

  void LoadMask(std::string);

  void StrokeUpdated(vtkPolyData* path, bool closed);
  void UpdateMaskFromStroke(vtkPolyData* path, bool closed);
  void UpdateColorPropagationLineFromStroke(vtkPolyData* polyDataPath);

  // A class to do the main computations in a separate thread so we can display a marquee progress bar.
  //CProgressThread ProgressThread;

  // Our scribble interactor style
  vtkSmartPointer<vtkScribbleInteractorStyle> ScribbleInteractorStyle;

  // The main class. This will be instantiated as a StructurePropagation after the user selects which type of image to open.
  StructurePropagation<TImage>* StructurePropagationFilter;

  // The input image actors
  vtkSmartPointer<vtkImageActor> OriginalImageActor;
  vtkSmartPointer<vtkImageActor> MaskImageActor;

  // The output image actors
  vtkSmartPointer<vtkImageActor> ResultActor;

  // The renderers
  vtkSmartPointer<vtkRenderer> LeftRenderer;
  vtkSmartPointer<vtkRenderer> RightRenderer;

  // Refresh both renderers and render windows
  void Refresh();

  // Allows the background color to be changed
  double BackgroundColor[3];

  // We set this when the image is opened. We sometimes need to know how big the image is.
  itk::ImageRegion<2> ImageRegion;

  // Data
  UnsignedCharScalarImageType::Pointer Mask;
  std::vector<itk::Index<2> > ColorPropagationLine;

  // Data, mapper, and actor for the selections
  vtkSmartPointer<vtkPolyData> ColorPropagationPathPolyData;
  vtkSmartPointer<vtkPolyDataMapper> ColorPropagationPathMapper;
  vtkSmartPointer<vtkActor> ColorPropagationPathActor;

};

#include "InnerWidget.txx"

#endif
