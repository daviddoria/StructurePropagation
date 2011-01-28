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

/*
 * This class is responsible for the user interaction with the input image.
 * A vtkImageTracerWidget does most of the work, but this class appends and maintains
 * the selections.
*/

#ifndef vtkScribbleInteractorStyle_H
#define vtkScribbleInteractorStyle_H

#include <vtkImageTracerWidget.h>
#include <vtkInteractorStyleImage.h> // superclass
#include <vtkSmartPointer.h>

#include "itkImageRegion.h"

#include "Types.h"

class vtkImageActor;
class vtkImageData;
class vtkPolyData;

class vtkScribbleInteractorStyle : public vtkInteractorStyleImage
{
public:
  static vtkScribbleInteractorStyle* New();
  vtkTypeMacro(vtkScribbleInteractorStyle, vtkInteractorStyleImage);

  vtkScribbleInteractorStyle();

  int GetStrokeType();
  enum SELECTION {COLOR};

  // A named function to set the SelectionType
  void SetInteractionModeToColor();

  // Get the strokes
  std::vector<itk::Index<2> > GetColorStrokes();

  // Clear all strokes
  void ClearStrokes();

  // Connect the tracer to the interactor, etc.
  void InitializeTracer(vtkImageActor* imageActor);

  // Get an image the size of the input image with a black background
  // and white pixels indicating the stroke
  void GetStrokeImage(UnsignedCharScalarImageType::Pointer image);

  // Get an image the size of the input image with a black background
  // and white pixels indicating a dilated representation of the stroke
  void GetDilatedStrokeImage(UnsignedCharScalarImageType::Pointer image);

private:
  void Refresh();

  // Update the selection when the EndInteraction event is fired.
  void CatchWidgetEvent(vtkObject* caller, long unsigned int eventId, void* callData);

  // The type of stroke currently selected.
  int StrokeType;

  // The widget which does most of the work.
  vtkSmartPointer<vtkImageTracerWidget> Tracer;

  // Keep track of the pixels the user selected.
  std::vector<itk::Index<2> > ColorStrokes;

  // Data, mapper, and actor for the selections
  vtkSmartPointer<vtkPolyData> ColorStrokePolyData;

  vtkSmartPointer<vtkPolyDataMapper> ColorStrokeMapper;

  vtkSmartPointer<vtkActor> ColorStrokeActor;

  // We want to be able to access images of the strokes
  itk::ImageRegion<2> ImageRegion;
};

#endif