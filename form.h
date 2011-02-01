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

/* This is the main GUI class of this project. It is a QMainWindow
 * so that we can use a File menu. It contains an instance of our main functional
 * class ImageGraphCutBase and our custom scribble interactor style vtkGraphCutInteractorStyle.
 * It also contains a CProgressThread so that we can display a progress bar in marquee
 * mode during long computations.
*/

#ifndef FILEMENUFORM_H
#define FILEMENUFORM_H

// Qt
#include "ui_form.h"

// ITK
#include "itkCovariantVector.h"
#include "itkImage.h"
#include "itkImageRegionConstIteratorWithIndex.h"

// Custom
#include "vtkScribbleInteractorStyle.h"
//#include "ImageGraphCutBase.h"
//#include "ImageGraphCut.h"
#include "ProgressThread.h"
#include "StructurePropagation.h"

// VTK
#include <vtkSmartPointer.h>
#include <vtkImageData.h>

// Forward declarations
class vtkImageActor;
class vtkRenderer;
class vtkImageData;

class Form : public QMainWindow, private Ui::MainWindow
{
Q_OBJECT
public:
  Form(QWidget *parent = 0);

public slots:
  // Menu items
  void actionOpen_Color_Image_triggered();
  void actionOpen_Grayscale_Image_triggered();
  void actionFlip_Image_triggered();
  void actionSave_Result_triggered();

  // Buttons, radio buttons, and sliders
  void btnClearStrokes_clicked();
  void btnSaveStrokes_clicked();
  void btnPropagate_clicked();
  void btnLoadMask_clicked();

  // These slots handle running the progress bar while the computations are done in a separate thread.
  void StartProgressSlot();
  void StopProgressSlot();

  // Testing
  void btnExtractPatches_clicked();

protected:

  // A class to do the main computations in a separate thread so we can display a marquee progress bar.
  CProgressThread ProgressThread;

  // Use a QFileDialog to get a filename, then open the specified file as a greyscale or color image, depending on which type the user has specified through the file menu.
  template<typename TImageType>
  void OpenFile();

  // Our scribble interactor style
  vtkSmartPointer<vtkScribbleInteractorStyle> ScribbleInteractorStyle;

  // The main class. This will be instantiated as a StructurePropagation after the user selects which type of image to open.
  StructurePropagationBase* StructurePropagationFilter;

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

};

#endif
