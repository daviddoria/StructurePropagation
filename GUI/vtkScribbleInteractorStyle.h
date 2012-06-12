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
 * This class emits a signal when the users scribbles on the image.
 */

#ifndef vtkScribbleInteractorStyle_H
#define vtkScribbleInteractorStyle_H

// VTK
#include <vtkImageTracerWidget.h>
#include <vtkInteractorStyleImage.h> // superclass
#include <vtkSmartPointer.h>

// Boost
#include <boost/signals2/signal.hpp>

class vtkImageActor;

class vtkScribbleInteractorStyle : public vtkInteractorStyleImage
{
public:
  // Constructors
  static vtkScribbleInteractorStyle* New();
  vtkTypeMacro(vtkScribbleInteractorStyle, vtkInteractorStyleImage);
  vtkScribbleInteractorStyle();

  // A signal to indicate that something has changed.
  boost::signals2::signal<void (vtkPolyData*, bool)> StrokeUpdated;

  // Connect the tracer to the interactor, etc.
  void InitializeTracer(vtkImageActor* imageActor);

private:

  void ClearTracer();

  // Update the selection when the EndInteraction event is fired.
  void CatchWidgetEvent(vtkObject* caller, long unsigned int eventId, void* callData);

  // The widget which does most of the work.
  vtkSmartPointer<vtkImageTracerWidget> Tracer;

};

#endif