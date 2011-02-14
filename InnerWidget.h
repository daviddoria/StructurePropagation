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
#include "ComputationThread.h"
#include "InnerWidgetObject.h"

// VTK
#include <vtkSmartPointer.h>

// Boost.Signals
#include <boost/signals2/signal.hpp>
#include <boost/bind.hpp>

#include <string>

// Forward declarations
class vtkImageActor;
class vtkRenderer;

template <typename TImage>
class InnerWidget : public InnerWidgetObject
{
public:
  InnerWidget(QWidget *parent = 0);

  // Functions
  void OpenFile();

  void chkFlip_clicked();

  // Buttons
  void btnClearStrokes_clicked();
  void btnSaveStrokes_clicked();
  void btnPropagate_clicked();
  void btnLoadMask_clicked();
  void btnSaveResult_clicked();
  void btnScreenshot_clicked();

  // These slots handle running the progress bar while the computations are done in a separate thread.
  void StartProgressSlot();
  void StopProgressSlot();

  // Check boxes
  void chkScale_clicked();
  void chkShowPaths_clicked();
  void chkShowResult_clicked();
  void chkShowOriginal_clicked();
  void chkShowMask_clicked();

  // Other GUI elements
  void sldPatchRadius_valueChanged();

protected:

  // Functions
  void ConnectSignalsAndSlots();
  void FlipImage();
  void SaveResult();
  void LoadMask(std::string filename);

  void DisplayTransparencyMaskedImage(typename TImage::Pointer image);
  void DisplayScaledImage(typename TImage::Pointer image);

  // Handle things associated with the propagation path
  void StrokeUpdated(vtkPolyData* path, bool closed);
  void UpdateMaskFromStroke(vtkPolyData* path, bool closed);
  void UpdateColorPropagationLineFromStroke(vtkPolyData* polyDataPath);

  // A class to do the main computations in a separate thread so we can display a marquee progress bar.
  ComputationThread<TImage> PropagateThread;

  // Scribble interactor style
  vtkSmartPointer<vtkScribbleInteractorStyle> ScribbleInteractorStyle;
  vtkSmartPointer<vtkImageActor> ScribbleCanvasActor;
  vtkSmartPointer<vtkImageData> ScribbleCanvas;
  void CreateScribbleCanvas();

  StructurePropagation<TImage> StructurePropagationFilter;

  // The input image actors
  vtkSmartPointer<vtkImageActor> OriginalImageActor;
  vtkSmartPointer<vtkImageActor> MaskImageActor;
  vtkSmartPointer<vtkImageActor> ResultActor;

  // The renderers
  vtkSmartPointer<vtkRenderer> Renderer;

  // Refresh both renderers and render windows
  void Refresh();
  bool NeverRendered;

  // Allows the background color to be changed
  double BackgroundColor[3];

  // Data
  UnsignedCharScalarImageType::Pointer Mask;
  //std::vector<itk::Index<2> > ColorPropagationLine;
  std::set<itk::Index<2> > ColorPropagationLine;
  std::set<itk::Index<2> > PropagationLineIntersections;
  itk::ImageRegion<2> ImageRegion; // This is set when the image is opened. We sometimes need to know how big the image is.

  // Objects for the propagation path.
  vtkSmartPointer<vtkPolyData> ColorPropagationPathPolyData;
  vtkSmartPointer<vtkPolyDataMapper> ColorPropagationPathMapper;
  vtkSmartPointer<vtkActor> ColorPropagationPathActor;

};

#include "InnerWidget.txx"

#endif
