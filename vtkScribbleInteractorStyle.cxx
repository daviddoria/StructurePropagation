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

#include "vtkScribbleInteractorStyle.h"

#include <vtkActor.h>
#include <vtkAppendPolyData.h>
#include <vtkCallbackCommand.h>
#include <vtkCommand.h>
#include <vtkImageActor.h>
#include <vtkImageData.h>
#include <vtkImageStencilToImage.h>
#include <vtkImageTracerWidget.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkProperty.h>
#include <vtkObjectFactory.h>
#include <vtkRendererCollection.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>

#include "Helpers.h"

vtkStandardNewMacro(vtkScribbleInteractorStyle);

vtkScribbleInteractorStyle::vtkScribbleInteractorStyle()
{
  // Initializations
  this->Tracer = vtkSmartPointer<vtkImageTracerWidget>::New();
  this->Tracer->GetLineProperty()->SetLineWidth(5);
  //this->Tracer->HandleMiddleMouseButtonOff();

  // Update the selection when the EndInteraction event is fired.
  this->Tracer->AddObserver(vtkCommand::EndInteractionEvent, this, &vtkScribbleInteractorStyle::CatchWidgetEvent);

}

void vtkScribbleInteractorStyle::InitializeTracer(vtkImageActor* imageActor)
{
  this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(imageActor);
  this->Tracer->SetInteractor(this->Interactor);
  this->Tracer->SetViewProp(imageActor);
  this->Tracer->ProjectToPlaneOn();

  this->Tracer->On();
}

void vtkScribbleInteractorStyle::CatchWidgetEvent(vtkObject* caller, long unsigned int eventId, void* callData)
{
  // Get the tracer object (this is the object that triggered this event)
  vtkImageTracerWidget* tracer =
    static_cast<vtkImageTracerWidget*>(caller);

  vtkSmartPointer<vtkPolyData> path =
    vtkSmartPointer<vtkPolyData>::New();
  tracer->GetPath(path);

  this->StrokeUpdated(path, tracer->IsClosed());

  ClearTracer();

};


void vtkScribbleInteractorStyle::ClearTracer()
{
  // "Clear" the tracer. We must rely on the foreground and background actors to maintain the appropriate colors.
  // If we did not clear the tracer, if we draw a foreground stroke (green) then switch to background mode, the last stoke would turn
  // red until we finished drawing the next stroke.
  vtkSmartPointer<vtkPoints> emptyPoints =
    vtkSmartPointer<vtkPoints>::New();
  emptyPoints->InsertNextPoint(0, 0, 0);
  emptyPoints->InsertNextPoint(0, 0, 0);

  this->Tracer->InitializeHandles(emptyPoints);
  this->Tracer->Modified();

  this->Refresh();
}

void vtkScribbleInteractorStyle::Refresh()
{
  this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->Render();
  this->Interactor->GetRenderWindow()->Render();
}
