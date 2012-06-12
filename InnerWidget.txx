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

// ITK
#include <itkCastImageFilter.h>
#include <itkCovariantVector.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkLineIterator.h>
#include <itkNthElementImageAdaptor.h>
#include "itkNumericTraits.h"

// VTK
#include <vtkAppendPolyData.h>
#include <vtkCamera.h>
#include <vtkImageActor.h>
#include <vtkImageData.h>
#include <vtkImageMapToColors.h>
#include <vtkImageStencilToImage.h>
#include <vtkInteractorStyleImage.h>
#include <vtkLookupTable.h>
#include <vtkPNGWriter.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkRenderWindow.h>
#include <vtkSmartPointer.h>
#include <vtkWindowToImageFilter.h>

// Qt
#include <QFileDialog>
#include <QLineEdit>
#include <QMessageBox>

// Custom
#include "Helpers.h"

template<typename TImage>
InnerWidget<TImage>::InnerWidget(QWidget *parent)
{
  this->NeverRendered = true;

  // Set the progress bar to marquee mode
  this->progressBar->setMinimum(0);
  this->progressBar->setMaximum(0);
  this->progressBar->hide();

  this->BackgroundColor[0] = 0;
  this->BackgroundColor[1] = 0;
  this->BackgroundColor[2] = .5;

  double initialCameraUp[3];
  initialCameraUp[0] = 0;
  initialCameraUp[1] = 1;
  initialCameraUp[2] = 0;

  // Instantiations
  this->OriginalImageActor = vtkSmartPointer<vtkImageActor>::New();
  this->MaskImageActor = vtkSmartPointer<vtkImageActor>::New();
  this->ResultActor = vtkSmartPointer<vtkImageActor>::New();

  this->ScribbleCanvasActor = vtkSmartPointer<vtkImageActor>::New();
  this->ScribbleCanvas = vtkSmartPointer<vtkImageData>::New();

  // Add renderers - we flip the image by changing the camera view up because of the conflicting conventions used by ITK and VTK
  this->Renderer = vtkSmartPointer<vtkRenderer>::New();
  this->Renderer->GradientBackgroundOn();
  this->Renderer->SetBackground(this->BackgroundColor);
  this->Renderer->SetBackground2(1,1,1);
  this->Renderer->GetActiveCamera()->SetViewUp(initialCameraUp);
  this->qvtkWidget->GetRenderWindow()->AddRenderer(this->Renderer);

  // Setup interactor style
  this->ScribbleInteractorStyle = vtkSmartPointer<vtkScribbleInteractorStyle>::New();
  this->qvtkWidget->GetInteractor()->SetInteractorStyle(this->ScribbleInteractorStyle);

  this->Mask = UnsignedCharScalarImageType::New();

  this->ColorPropagationPathPolyData = vtkSmartPointer<vtkPolyData>::New();
  this->ColorPropagationPathMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  this->ColorPropagationPathMapper->SetInputConnection(this->ColorPropagationPathPolyData->GetProducerPort());
  this->ColorPropagationPathActor = vtkSmartPointer<vtkActor>::New();
  this->ColorPropagationPathActor->SetMapper(this->ColorPropagationPathMapper);
  this->ColorPropagationPathActor->GetProperty()->SetLineWidth(4);
  this->ColorPropagationPathActor->GetProperty()->SetColor(0,1,0);
  //this->Renderer->AddActor(this->ColorPropagationPathActor);

  ConnectSignalsAndSlots();
}

template <typename TImage>
void InnerWidget<TImage>::ConnectSignalsAndSlots()
{
  this->ScribbleInteractorStyle->StrokeUpdated.connect(boost::bind(&InnerWidget::StrokeUpdated, this, _1, _2));

  connect(&PropagateThread, SIGNAL(StartProgressSignal()), this, SLOT(StartProgressSlot()), Qt::QueuedConnection);
  connect(&PropagateThread, SIGNAL(StopProgressSignal()), this, SLOT(StopProgressSlot()), Qt::QueuedConnection);

  connect(this->chkScale, SIGNAL(clicked()), this, SLOT(chkScale_clicked()));
  connect(this->chkFlip, SIGNAL(clicked()), this, SLOT(chkFlip_clicked()));
  connect(this->chkShowPaths, SIGNAL(clicked()), this, SLOT(chkShowPaths_clicked()));
  connect(this->chkShowResult, SIGNAL(clicked()), this, SLOT(chkShowResult_clicked()));
  connect(this->chkShowOriginal, SIGNAL(clicked()), this, SLOT(chkShowOriginal_clicked()));
  connect(this->chkShowMask, SIGNAL(clicked()), this, SLOT(chkShowMask_clicked()));

  connect(this->btnSaveResult, SIGNAL(clicked()), this, SLOT(btnSaveResult_clicked()));
  connect(this->btnScreenshot, SIGNAL(clicked()), this, SLOT(btnScreenshot_clicked()));

  connect( this->btnPropagate, SIGNAL( clicked() ), this, SLOT(btnPropagate_clicked()));
  connect( this->btnLoadMask, SIGNAL( clicked() ), this, SLOT(btnLoadMask_clicked()));
  connect( this->btnClearStrokes, SIGNAL( clicked() ), this, SLOT(btnClearStrokes_clicked()));
  connect( this->btnSaveStrokes, SIGNAL( clicked() ), this, SLOT(btnSaveStrokes_clicked()));
}

template <typename TImage>
void InnerWidget<TImage>::btnScreenshot_clicked()
{
  vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter =
    vtkSmartPointer<vtkWindowToImageFilter>::New();
  //windowToImageFilter->SetInput(renderWindow);
  windowToImageFilter->SetInput(this->qvtkWidget->GetInteractor()->GetRenderWindow());
  //windowToImageFilter->SetMagnification(3); //set the resolution of the output image (3 times the current resolution of vtk render window)
  //windowToImageFilter->SetInputBufferTypeToRGBA(); //also record the alpha (transparency) channel
  windowToImageFilter->Update();

  vtkSmartPointer<vtkPNGWriter> writer =
    vtkSmartPointer<vtkPNGWriter>::New();
  writer->SetFileName("screenshot.png");
  writer->SetInput(windowToImageFilter->GetOutput());
  writer->Write();
}

template <typename TImage>
void InnerWidget<TImage>::chkShowPaths_clicked()
{
  //this->ColorPropagationPathActor->SetVisibility(this->chkShowPaths->isChecked());
  this->Refresh();
  std::cout << "Set path visibility to " << this->chkShowPaths->isChecked() << std::endl;
}

template <typename TImage>
void InnerWidget<TImage>::chkShowResult_clicked()
{
  //this->ResultActor->SetVisibility(this->chkShowResult->isChecked());
  this->Refresh();
  std::cout << "Set result visibility to " << this->chkShowResult->isChecked() << std::endl;
}

template <typename TImage>
void InnerWidget<TImage>::chkShowOriginal_clicked()
{
  //this->OriginalImageActor->SetVisibility(this->chkShowOriginal->isChecked());
  this->Refresh();
  std::cout << "Set original image visibility to " << this->chkShowOriginal->isChecked() << std::endl;
}

template <typename TImage>
void InnerWidget<TImage>::chkShowMask_clicked()
{
  //this->MaskImageActor->SetVisibility(this->chkShowMask->isChecked());
  this->Refresh();
  std::cout << "Set mask visibility to " << this->chkShowMask->isChecked() << std::endl;
}


template <typename TImage>
void InnerWidget<TImage>::chkScale_clicked()
{
  /*
  if(this->chkScale->isChecked())
    {
    DisplayScaledImage(this->Image);
    }
  else
    {
    DisplayImage(this->Image);
    }
  */
}

template <typename TImage>
void InnerWidget<TImage>::StartProgressSlot()
{
  // Connected to the StartProgressSignal of the ProgressThread member
  std::cout << "StartProgressSlot: Called!" << std::endl;
  this->progressBar->show();
  QCursor waitCursor;
  waitCursor.setShape(Qt::WaitCursor);
  this->setCursor(waitCursor);
}

template <typename TImage>
void InnerWidget<TImage>::StopProgressSlot()
{
  // Display result

  //std::cout << "StopProgressSlot: Done!" << std::endl;
  this->progressBar->hide();
  QCursor arrowCursor;
  arrowCursor.setShape(Qt::ArrowCursor);
  this->setCursor(arrowCursor);

  //DisplayImage(this->StructurePropagationFilter.GetOutputImage());
  DisplayTransparencyMaskedImage(this->StructurePropagationFilter.GetOutputImage());

  this->SaveResult();
}

template<typename TImage>
void InnerWidget<TImage>::SaveResult()
{
  QFileInfo myFile(this->txtResult->text().toStdString().c_str());
  if(myFile.suffix().toStdString().compare(".png"))
    {
    Helpers::CastAndWriteImage<TImage>(this->StructurePropagationFilter.GetOutputImage(), this->txtResult->text().toStdString());
    }
  else
    {
    Helpers::WriteImage<TImage>(this->StructurePropagationFilter.GetOutputImage(), this->txtResult->text().toStdString());
    }
}

template<typename TImage>
void InnerWidget<TImage>::DisplayScaledImage(typename TImage::Pointer image)
{
  // Maybe this could be replaced with a "Magnitude of vector image" filter" (so it makes more sense with TImage more than 1 component) - an implicit "vector -> grayscale" computation - AbsImageFilter perhaps?
  typedef itk::NthElementImageAdaptor<TImage, float> ImageAdaptorType;
  typename ImageAdaptorType::Pointer adaptor = ImageAdaptorType::New();
  adaptor->SelectNthElement(0);
  adaptor->SetImage(image);

  typedef itk::RescaleIntensityImageFilter< ImageAdaptorType, UnsignedCharScalarImageType > RescaleFilterType;
  typename RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  rescaleFilter->SetInput(adaptor);
  rescaleFilter->SetOutputMinimum(0);
  rescaleFilter->SetOutputMaximum(255);
  rescaleFilter->Update();

  // Convert the ITK image to a VTK image and display it
  vtkSmartPointer<vtkImageData> VTKImage =
    vtkSmartPointer<vtkImageData>::New();
  Helpers::ITKImageToVTKImage<UnsignedCharScalarImageType>(rescaleFilter->GetOutput(), VTKImage);

  this->OriginalImageActor->SetInput(VTKImage);

}


template<typename TImage>
void InnerWidget<TImage>::DisplayTransparencyMaskedImage(typename TImage::Pointer image)
{
  // Convert the ITK image to a VTK image and display it
  vtkSmartPointer<vtkImageData> VTKImage =
    vtkSmartPointer<vtkImageData>::New();
  Helpers::ITKImageToVTKImage<TImage>(image, VTKImage);

  /*
  vtkSmartPointer<vtkImageData> transparencyMaskedImage =
    vtkSmartPointer<vtkImageData>::New();
  Helpers::ApplyTransparencyMask(VTKImage, this->Mask, transparencyMaskedImage);

  this->ResultActor->SetInput(transparencyMaskedImage);
  */

  this->ResultActor->SetInput(VTKImage);

  this->Renderer->AddActor(this->ResultActor);
  this->Refresh();
}

template<typename TImage>
void InnerWidget<TImage>::LoadMask(std::string filename)
{

  // Read file
  itk::ImageFileReader<UnsignedCharScalarImageType>::Pointer reader =
    itk::ImageFileReader<UnsignedCharScalarImageType>::New();
  reader->SetFileName(filename);
  reader->Update();

  this->Mask->Graft(reader->GetOutput());

  vtkSmartPointer<vtkImageData> VTKMaskImage =
    vtkSmartPointer<vtkImageData>::New();
  Helpers::ITKImageToVTKImage<UnsignedCharScalarImageType>(reader->GetOutput(), VTKMaskImage);

  vtkSmartPointer<vtkLookupTable> lookupTable =
    vtkSmartPointer<vtkLookupTable>::New();
  lookupTable->SetNumberOfTableValues(2);
  lookupTable->SetRange(0.0,1.0);
  lookupTable->SetTableValue( 0, 0.0, 0.0, 0.0, 0.0 ); //label 0 is transparent
  lookupTable->SetTableValue( 1, 1.0, 1.0, 1.0, 1.0 ); //label 1 is opaque and white
  lookupTable->Build();

  vtkSmartPointer<vtkImageMapToColors> mapTransparency =
    vtkSmartPointer<vtkImageMapToColors>::New();
  mapTransparency->SetLookupTable(lookupTable);
  mapTransparency->SetInput(VTKMaskImage);
  mapTransparency->PassAlphaToOutputOn();

  this->MaskImageActor->SetInput(mapTransparency->GetOutput());

  //this->qvtkWidget->GetInteractor()->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(this->MaskImageActor);

}

template<typename TImage>
void InnerWidget<TImage>::btnLoadMask_clicked()
{
  // Get a filename to open
  QString filename = QFileDialog::getOpenFileName(this, "Open Mask", ".", "Image Files (*.png *.bmp)");

  if(filename.isEmpty())
    {
    std::cerr << "No file selected!" << std::endl;
    return;
    }

  LoadMask(filename.toStdString());
}

template<typename TImage>
void InnerWidget<TImage>::btnClearStrokes_clicked()
{
  this->ColorPropagationLine.clear();
  this->ColorPropagationPathPolyData->Reset();
  this->ColorPropagationPathPolyData->Squeeze();
  this->ColorPropagationPathMapper->Update();
  this->Refresh();
}

template<typename TImage>
void InnerWidget<TImage>::btnSaveStrokes_clicked()
{
  QString directoryName = QFileDialog::getExistingDirectory(this,
    "Chose a directory to save the stroke images", ".", QFileDialog::ShowDirsOnly);

  std::string colorStrokesFilename = QDir(directoryName).absoluteFilePath("ColorStrokes.png").toStdString();

  std::cout << "Writing to " << colorStrokesFilename << std::endl;

  // Get image of stroke from scribble interactor style
  UnsignedCharScalarImageType::Pointer colorStrokeImage = UnsignedCharScalarImageType::New();
  Helpers::IndicesToBinaryImage(this->ColorPropagationLine, this->Mask->GetLargestPossibleRegion(), colorStrokeImage);

  Helpers::WriteImage<UnsignedCharScalarImageType>(colorStrokeImage, colorStrokesFilename);

}

template<typename TImage>
void InnerWidget<TImage>::chkFlip_clicked()
{
  this->FlipImage();
}

template<typename TImage>
void InnerWidget<TImage>::btnPropagate_clicked()
{
  std::cout << "btnPropagate_clicked: Number of points: " << this->ColorPropagationLine.size() << std::endl;

  this->StructurePropagationFilter.ClearEverything();

  this->StructurePropagationFilter.SetPatchRadius(this->scrPatchRadius->value());
  this->StructurePropagationFilter.SetSourceLineWidth(this->scrSourceLineWidth->value());
  this->StructurePropagationFilter.SetMask(this->Mask);
  this->StructurePropagationFilter.SetPropagationLine(this->ColorPropagationLine);
  this->StructurePropagationFilter.SetPropagationLineIntersections(this->PropagationLineIntersections);

  // Perform the propagation in a separate thread
  this->PropagateThread.Propagation = &(this->StructurePropagationFilter);
  this->PropagateThread.start();
}


template<typename TImage>
void InnerWidget<TImage>::FlipImage()
{
  double up[3];
  this->Renderer->GetActiveCamera()->GetViewUp(up);
  up[1] *= -1;
  this->Renderer->GetActiveCamera()->SetViewUp(up);

  double pos[3];
  this->Renderer->GetActiveCamera()->GetPosition(pos);
  pos[2] *= -1;
  this->Renderer->GetActiveCamera()->SetPosition(pos);

  this->Refresh();

}


template<typename TImage>
void InnerWidget<TImage>::btnSaveResult_clicked()
{
  this->SaveResult();
}


template <typename TImage>
void InnerWidget<TImage>::OpenFile()
{
  // Use a QFileDialog to get a filename, then open the specified file as a greyscale or color image, depending on which type the user has specified through the file menu.

  // Get a filename to open
  QString filename = QFileDialog::getOpenFileName(this,
     //tr("Open Image"), "/media/portable/Projects/src/InteractiveImageGraphCutSegmentation/data", tr("Image Files (*.png *.bmp *.mhd)"));
      tr("Open Image"), ".", tr("Image Files (*.png *.bmp *.mhd)"));

  if(filename.isEmpty())
    {
    return;
    }

  // Read file
  typename itk::ImageFileReader<TImage>::Pointer reader = itk::ImageFileReader<TImage>::New();

  reader->SetFileName(filename.toStdString());
  reader->Update();

  this->ImageRegion = reader->GetOutput()->GetLargestPossibleRegion();

  if(this->chkScale->isChecked())
    {
    //DisplayScaledImage(reader->GetOutput());
    }
  else
    {
    vtkSmartPointer<vtkImageData> VTKImage =
      vtkSmartPointer<vtkImageData>::New();
    Helpers::ITKImageToVTKImage<TImage>(reader->GetOutput(), VTKImage);
    this->OriginalImageActor->SetInput(VTKImage);
    }

  this->StructurePropagationFilter.SetImage(reader->GetOutput());

  // Create mask filename from image filename
  QFileInfo maskFileInfo(filename);
  std::stringstream ssMaskFile;
  ssMaskFile << maskFileInfo.absolutePath().toStdString() << "/" << maskFileInfo.baseName().toStdString() << "Mask.png";
  std::cout << "Constructed default mask filename as " << ssMaskFile.str() << std::endl;

  // If the mask file exists, load it.
  QFile myFile(ssMaskFile.str().c_str());
  if(myFile.exists())
    {
    std::cout << "Loading mask " << ssMaskFile.str() << std::endl;
    LoadMask(ssMaskFile.str());
    }
  else
    {
    std::cout << "Mask " << ssMaskFile.str() << " not found!" << std::endl;
    }

  CreateScribbleCanvas();

  this->NeverRendered = true;

  this->Refresh();

  // Optionally flip the image (if you are frequently using images which must be flipped
  // FlipImage();

}

template<typename TImage>
void InnerWidget<TImage>::Refresh()
{
  //std::cout << "Refresh()" << std::endl;

  // Remove all actors
  this->Renderer->RemoveActor(ColorPropagationPathActor);
  this->Renderer->RemoveActor(MaskImageActor);
  this->Renderer->RemoveActor(OriginalImageActor);
  this->Renderer->RemoveActor(ResultActor);
  this->Renderer->RemoveActor(ScribbleCanvasActor);

  // Add all actors back if they are enabled (order is important here)
  if(this->chkShowOriginal->isChecked())
    {
    this->Renderer->AddActor(OriginalImageActor);
    }
  if(this->chkShowMask->isChecked())
    {
    this->Renderer->AddActor(MaskImageActor);
    }
  if(this->chkShowResult->isChecked())
    {
    this->Renderer->AddActor(ResultActor);
    }
  if(this->chkShowPaths->isChecked())
    {
    this->Renderer->AddActor(ColorPropagationPathActor);
    }

  this->Renderer->AddActor(ScribbleCanvasActor);

  if(this->NeverRendered)
    {
    this->Renderer->ResetCamera();
    this->NeverRendered = false;
    }

  this->Renderer->Render();
  this->qvtkWidget->GetRenderWindow()->Render();
  this->qvtkWidget->GetInteractor()->Render();

  this->ScribbleInteractorStyle->InitializeTracer(this->ScribbleCanvasActor);
}

template<typename TImage>
void InnerWidget<TImage>::CreateScribbleCanvas()
{
  this->ScribbleCanvas->SetNumberOfScalarComponents(4);
  this->ScribbleCanvas->SetScalarTypeToUnsignedChar();
  int dims[3];
  this->OriginalImageActor->GetInput()->GetDimensions(dims);
  this->ScribbleCanvas->SetDimensions(dims);
  this->ScribbleCanvas->AllocateScalars();

  for (int y = 0; y < dims[1]; y++)
    {
    for (int x = 0; x < dims[0]; x++)
      {
      unsigned char* pixel = static_cast<unsigned char*>(this->ScribbleCanvas->GetScalarPointer(x,y,0));
      pixel[3] = 0; // transparent
      }
    }

  this->ScribbleCanvasActor->SetInput(this->ScribbleCanvas);

  //this->ScribbleInteractorStyle->InitializeTracer(this->ScribbleCanvasActor);
}

template<typename TImage>
void InnerWidget<TImage>::StrokeUpdated(vtkPolyData* path, bool closed)
{
  UpdateColorPropagationLineFromStroke(path);
}

template<typename TImage>
void InnerWidget<TImage>::UpdateColorPropagationLineFromStroke(vtkPolyData* polyDataPath)
{
  //std::cout << "Line size before update: " << this->ColorPropagationLine.size() << std::endl;
  //std::cout << "Intersection size before update: " << this->PropagationLineIntersections.size() << std::endl;

  vtkSmartPointer<vtkAppendPolyData> appendFilter =
    vtkSmartPointer<vtkAppendPolyData>::New();
  appendFilter->AddInputConnection(this->ColorPropagationPathPolyData->GetProducerPort());
  appendFilter->AddInputConnection(polyDataPath->GetProducerPort());
  appendFilter->Update();
  this->ColorPropagationPathPolyData->ShallowCopy(appendFilter->GetOutput());

  std::vector<itk::Index<2> > path = Helpers::PolyDataToPixelList(polyDataPath);
  typedef std::pair<std::set<itk::Index<2> >::iterator, bool> ReturnType;
  for(unsigned int i = 0; i < path.size(); i++)
    {
    //std::cout << "path pixel " << i << " : " << path[i] << std::endl;

    ReturnType inserted = this->ColorPropagationLine.insert(path[i]);
    if(inserted.second == false)
      {
      this->PropagationLineIntersections.insert(path[i]);
      //std::cout << "UpdateColorPropagationLineFromStroke: Intersection at " << path[i] << std::endl;
      }
    }

  //std::cout << "Size after update: " << this->ColorPropagationLine.size() << std::endl;
  std::cout << "Total intersections:: " << this->PropagationLineIntersections.size() << std::endl;
}

template<typename TImage>
void InnerWidget<TImage>::UpdateMaskFromStroke(vtkPolyData* path, bool closed)
{
  // We need a VTK image to get the spacing, origin, and extent only
  vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();
  Helpers::ITKImageToVTKImage<UnsignedCharScalarImageType>(this->Mask, image);

  vtkSmartPointer<vtkPolyDataToImageStencil> polyDataToImageStencil =
    vtkSmartPointer<vtkPolyDataToImageStencil>::New();
  polyDataToImageStencil->SetTolerance(0);
  polyDataToImageStencil->SetInputConnection(path->GetProducerPort());
  polyDataToImageStencil->SetOutputOrigin(image->GetOrigin());
  polyDataToImageStencil->SetOutputSpacing(image->GetSpacing());
  polyDataToImageStencil->SetOutputWholeExtent(image->GetExtent());
  polyDataToImageStencil->Update();

  vtkSmartPointer<vtkImageStencilToImage> imageStencilToImage =
    vtkSmartPointer<vtkImageStencilToImage>::New();
  imageStencilToImage->SetInputConnection(polyDataToImageStencil->GetOutputPort());
  imageStencilToImage->SetInsideValue(255);
  imageStencilToImage->Update();

  UnsignedCharScalarImageType::Pointer mask = UnsignedCharScalarImageType::New();
  Helpers::VTKImageToITKImage(imageStencilToImage->GetOutput(), mask);

}