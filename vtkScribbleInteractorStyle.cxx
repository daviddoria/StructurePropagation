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
#include <vtkImageTracerWidget.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkObjectFactory.h>
#include <vtkRendererCollection.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>

#include <itkImageRegionIterator.h>
#include <itkBresenhamLine.h>
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"

#include "Helpers.h"

vtkStandardNewMacro(vtkScribbleInteractorStyle);

vtkScribbleInteractorStyle::vtkScribbleInteractorStyle()
{
  // Initializations
  this->Tracer = vtkSmartPointer<vtkImageTracerWidget>::New();
  this->Tracer->GetLineProperty()->SetLineWidth(5);
  this->Tracer->HandleMiddleMouseButtonOff();

  // Foreground
  this->ColorStrokePolyData = vtkSmartPointer<vtkPolyData>::New();
  this->ColorStrokeMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  this->ColorStrokeActor = vtkSmartPointer<vtkActor>::New();
  this->ColorStrokeActor->SetMapper(this->ColorStrokeMapper);
  this->ColorStrokeActor->GetProperty()->SetLineWidth(4);
  this->ColorStrokeActor->GetProperty()->SetColor(0,1,0);
  this->ColorStrokeMapper->SetInputConnection(this->ColorStrokePolyData->GetProducerPort());

  // Update the selection when the EndInteraction event is fired.
  this->Tracer->AddObserver(vtkCommand::EndInteractionEvent, this, &vtkScribbleInteractorStyle::CatchWidgetEvent);

  // Defaults
  this->StrokeType = COLOR;
}

std::vector<itk::Index<2> > vtkScribbleInteractorStyle::GetColorStrokes()
{
  return this->ColorStrokes;
}

int vtkScribbleInteractorStyle::GetStrokeType()
{
  return this->StrokeType;
}

void vtkScribbleInteractorStyle::InitializeTracer(vtkImageActor* imageActor)
{
  this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(imageActor);
  this->Tracer->SetInteractor(this->Interactor);
  this->Tracer->SetViewProp(imageActor);
  this->Tracer->ProjectToPlaneOn();

  this->Tracer->On();

  itk::Index<2> start;
  start.Fill(0);
  itk::Size<2> size;
  int dims[3];
  imageActor->GetInput()->GetDimensions(dims);
  size[0] = dims[0];
  size[1] = dims[1];
  this->ImageRegion.SetSize(size);
  this->ImageRegion.SetIndex(start);
}

void vtkScribbleInteractorStyle::SetInteractionModeToColor()
{
  this->Tracer->GetLineProperty()->SetColor(0,1,0);
  this->StrokeType = COLOR;
}


void vtkScribbleInteractorStyle::CatchWidgetEvent(vtkObject* caller, long unsigned int eventId, void* callData)
{
  // Get the path from the tracer and append it to the appropriate selection

  this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(ColorStrokeActor);
  //this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(ForegroundSelectionActor);

  // Get the tracer object (this is the object that triggered this event)
  vtkImageTracerWidget* tracer =
    static_cast<vtkImageTracerWidget*>(caller);

  // Get the points in the selection
  vtkSmartPointer<vtkPolyData> path =
    vtkSmartPointer<vtkPolyData>::New();
  tracer->GetPath(path);

  // Create a filter which will be used to combine the most recent selection with previous selections
  vtkSmartPointer<vtkAppendPolyData> appendFilter =
    vtkSmartPointer<vtkAppendPolyData>::New();
  appendFilter->AddInputConnection(path->GetProducerPort());

  std::vector<itk::Index<2> > newPoints = PolyDataToPixelList(path);
  //std::cout << newPoints.size() << " new points." << std::endl;

  // If we are in foreground mode, add the current selection to the foreground. Else, add it to the background.
  if(this->StrokeType == vtkScribbleInteractorStyle::COLOR)
    {
    appendFilter->AddInputConnection(this->ColorStrokePolyData->GetProducerPort());
    appendFilter->Update();
    this->ColorStrokePolyData->ShallowCopy(appendFilter->GetOutput());

    this->ColorStrokes.insert(this->ColorStrokes.end(), newPoints.begin(), newPoints.end());
    }

  //std::cout << this->ForegroundSelection.size() << " foreground poitns." << std::endl;
  //std::cout << this->BackgroundSelection.size() << " background poitns." << std::endl;

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

};

void vtkScribbleInteractorStyle::Refresh()
{
  this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->Render();
  this->Interactor->GetRenderWindow()->Render();
}

void vtkScribbleInteractorStyle::ClearStrokes()
{
  /*
   // I thought this would work...
  this->BackgroundSelection->Reset();
  this->BackgroundSelection->Squeeze();
  this->BackgroundSelection->Modified();

  this->ForegroundSelection->Reset();
  this->ForegroundSelection->Squeeze();
  this->ForegroundSelection->Modified();
  */

  // This seems like a silly way of emptying the polydatas...
  vtkSmartPointer<vtkPolyData> empytPolyData =
    vtkSmartPointer<vtkPolyData>::New();
  this->ColorStrokePolyData->ShallowCopy(empytPolyData);

  this->ColorStrokes.clear();

  this->Refresh();
}



void vtkScribbleInteractorStyle::GetStrokeImage(UnsignedCharScalarImageType::Pointer image)
{
  image->SetRegions(this->ImageRegion);
  image->Allocate();
  image->FillBuffer(itk::NumericTraits<UnsignedCharScalarImageType::PixelType>::Zero);

  IndicesToBinaryImage(this->ColorStrokes, image);
}


void vtkScribbleInteractorStyle::GetDilatedStrokeImage(UnsignedCharScalarImageType::Pointer image)
{
  image->SetRegions(this->ImageRegion);
  image->Allocate();
  image->FillBuffer(itk::NumericTraits<UnsignedCharScalarImageType::PixelType>::Zero);

  IndicesToBinaryImage(this->ColorStrokes, image);

  typedef itk::BinaryBallStructuringElement<UnsignedCharScalarImageType::PixelType,2> StructuringElementType;
  StructuringElementType structuringElement;
  structuringElement.SetRadius(2);
  structuringElement.CreateStructuringElement();

  typedef itk::BinaryDilateImageFilter <UnsignedCharScalarImageType, UnsignedCharScalarImageType, StructuringElementType>
          BinaryDilateImageFilterType;

  BinaryDilateImageFilterType::Pointer dilateFilter = BinaryDilateImageFilterType::New();
  dilateFilter->SetInput(image);
  dilateFilter->SetKernel(structuringElement);
  dilateFilter->Update();

  image->Graft(dilateFilter->GetOutput());
}