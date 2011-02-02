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

/*
This implementation is based on
Image Completion With Structure Propagation

*/

#include "form.h"
#include "Helpers.h"

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
#include <vtkCamera.h>
#include <vtkImageActor.h>
#include <vtkImageData.h>
#include <vtkImageMapToColors.h>
#include <vtkImageStencilToImage.h>
#include <vtkInteractorStyleImage.h>
#include <vtkLookupTable.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkRenderWindow.h>
#include <vtkSmartPointer.h>

// Qt
#include <QFileDialog>
#include <QLineEdit>
#include <QMessageBox>


Form::Form(QWidget *parent)
{
  // Setup the GUI and connect all of the signals and slots
  setupUi(this);
  connect( this->actionOpen_Color_Image, SIGNAL( triggered() ), this, SLOT(actionOpen_Color_Image_triggered()) );
  connect( this->actionOpen_Grayscale_Image, SIGNAL( triggered() ), this, SLOT(actionOpen_Grayscale_Image_triggered()) );

  // Buttons
  connect( this->btnPropagate, SIGNAL( clicked() ), this, SLOT(btnPropagate_clicked()));
  connect( this->btnLoadMask, SIGNAL( clicked() ), this, SLOT(btnLoadMask_clicked()));
  connect( this->btnClearStrokes, SIGNAL( clicked() ), this, SLOT(btnClearStrokes_clicked()));
  connect( this->btnSaveStrokes, SIGNAL( clicked() ), this, SLOT(btnSaveStrokes_clicked()));

  // Radio buttons
  //connect( this->radDrawHole, SIGNAL( clicked() ), this, SLOT(radDrawHole_clicked()) );
  //connect( this->radDrawPropagationLine, SIGNAL( clicked() ), this, SLOT(radDrawPropagationLine_clicked()) );

  // Menu
  connect( this->actionFlip_Image, SIGNAL( triggered() ), this, SLOT(actionFlip_Image_triggered()));
  connect( this->actionSave_Result, SIGNAL( triggered() ), this, SLOT(actionSave_Result_triggered()));

  // Thread
  connect(&ProgressThread, SIGNAL(StartProgressSignal()), this, SLOT(StartProgressSlot()), Qt::QueuedConnection);
  connect(&ProgressThread, SIGNAL(StopProgressSignal()), this, SLOT(StopProgressSlot()), Qt::QueuedConnection);

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

  // Add renderers - we flip the image by changing the camera view up because of the conflicting conventions used by ITK and VTK
  this->LeftRenderer = vtkSmartPointer<vtkRenderer>::New();
  this->LeftRenderer->GradientBackgroundOn();
  this->LeftRenderer->SetBackground(this->BackgroundColor);
  this->LeftRenderer->SetBackground2(1,1,1);
  this->LeftRenderer->GetActiveCamera()->SetViewUp(initialCameraUp);
  this->qvtkWidgetLeft->GetRenderWindow()->AddRenderer(this->LeftRenderer);

  this->RightRenderer = vtkSmartPointer<vtkRenderer>::New();
  this->RightRenderer->GradientBackgroundOn();
  this->RightRenderer->SetBackground(this->BackgroundColor);
  this->RightRenderer->SetBackground2(1,1,1);
  this->RightRenderer->GetActiveCamera()->SetViewUp(initialCameraUp);
  this->qvtkWidgetRight->GetRenderWindow()->AddRenderer(this->RightRenderer);

  // Setup right interactor style
  vtkSmartPointer<vtkInteractorStyleImage> interactorStyleImage =
    vtkSmartPointer<vtkInteractorStyleImage>::New();
  this->qvtkWidgetRight->GetInteractor()->SetInteractorStyle(interactorStyleImage);

  // Setup left interactor style
  this->ScribbleInteractorStyle = vtkSmartPointer<vtkScribbleInteractorStyle>::New();
  this->qvtkWidgetLeft->GetInteractor()->SetInteractorStyle(this->ScribbleInteractorStyle);

  this->StructurePropagationFilter = NULL;
  this->Mask = UnsignedCharScalarImageType::New();

  this->ColorPropagationPathPolyData = vtkSmartPointer<vtkPolyData>::New();
  this->ColorPropagationPathMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  this->ColorPropagationPathMapper->SetInputConnection(this->ColorPropagationPathPolyData->GetProducerPort());
  this->ColorPropagationPathActor = vtkSmartPointer<vtkActor>::New();
  this->ColorPropagationPathActor->SetMapper(this->ColorPropagationPathMapper);
  this->ColorPropagationPathActor->GetProperty()->SetLineWidth(4);
  this->ColorPropagationPathActor->GetProperty()->SetColor(0,1,0);
  this->LeftRenderer->AddActor(this->ColorPropagationPathActor);

  this->ScribbleInteractorStyle->StrokeUpdated.connect(boost::bind(&Form::StrokeUpdated, this, _1, _2));
}

void Form::btnLoadMask_clicked()
{
  if(!this->StructurePropagationFilter)
    {
    std::cerr << "Cannot load a mask until an image is loaded!" << std::endl;
    return;
    }

  // Get a filename to open
  QString filename = QFileDialog::getOpenFileName(this,
      tr("Open Mask"), ".", tr("Image Files (*.png *.bmp)"));

  if(filename.isEmpty())
    {
    std::cerr << "No file selected!" << std::endl;
    return;
    }

  // Read file
  itk::ImageFileReader<UnsignedCharScalarImageType>::Pointer reader =
    itk::ImageFileReader<UnsignedCharScalarImageType>::New();
  reader->SetFileName(filename.toStdString());
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

  //this->MaskImageActor->SetInput(VTKMaskImage);
  this->MaskImageActor->SetInput(mapTransparency->GetOutput());

  this->qvtkWidgetLeft->GetInteractor()->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(this->MaskImageActor);

  this->ScribbleInteractorStyle->InitializeTracer(this->MaskImageActor);
}

void Form::StartProgressSlot()
{
  // Connected to the StartProgressSignal of the ProgressThread member
  this->progressBar->show();
}

// Display result image with transparent background pixels
void Form::StopProgressSlot()
{
  // When the ProgressThread emits the StopProgressSignal, we need to display the result of the segmentation
#if 0
  // Convert the segmentation mask to a binary VTK image
  vtkSmartPointer<vtkImageData> VTKSegmentMask =
    vtkSmartPointer<vtkImageData>::New();
  ITKImagetoVTKImage<MaskImageType>(static_cast<ImageGraphCut<GrayscaleImageType>* >(this->GraphCut)->GetSegmentMask(), VTKSegmentMask);

  // Convert the image into a VTK image for display
  vtkSmartPointer<vtkImageData> VTKImage =
    vtkSmartPointer<vtkImageData>::New();
  if(this->GraphCut->GetPixelDimensionality() == 1)
    {
    ITKImagetoVTKImage<GrayscaleImageType>(static_cast<ImageGraphCut<GrayscaleImageType>* >(this->GraphCut)->GetMaskedOutput(), VTKImage);
    }
  else if(this->GraphCut->GetPixelDimensionality() == 3)
    {
    ITKImagetoVTKImage<ColorImageType>(static_cast<ImageGraphCut<ColorImageType>* >(this->GraphCut)->GetMaskedOutput(), VTKImage);
    }
  else
    {
    std::cerr << "This type of image (" << this->GraphCut->GetPixelDimensionality() << ") cannot be displayed!" << std::endl;
    exit(-1);
    }

  vtkSmartPointer<vtkImageData> VTKMaskedImage =
    vtkSmartPointer<vtkImageData>::New();
  MaskImage(VTKImage, VTKSegmentMask, VTKMaskedImage);

  // Remove the old output, set the new output and refresh everything
  this->ResultActor->SetInput(VTKMaskedImage);
  this->RightRenderer->RemoveAllViewProps();
  this->RightRenderer->AddActor(ResultActor);
  this->RightRenderer->ResetCamera();
  this->Refresh();

  this->progressBar->hide();
#endif
}


void Form::btnClearStrokes_clicked()
{
  this->ColorPropagationLine.clear();
}

void Form::btnSaveStrokes_clicked()
{
  QString directoryName = QFileDialog::getExistingDirectory(this,
    "Chose a directory to save the stroke images", ".", QFileDialog::ShowDirsOnly);

  std::string colorStrokesFilename = QDir(directoryName).absoluteFilePath("ColorStrokes.png").toStdString();

  std::cout << "Writing to " << colorStrokesFilename << std::endl;

  // Get image of stroke from scribble interactor style
  UnsignedCharScalarImageType::Pointer colorStrokeImage = UnsignedCharScalarImageType::New();
  Helpers::IndicesToBinaryImage(this->ColorPropagationLine, colorStrokeImage);

  typedef  itk::ImageFileWriter< UnsignedCharScalarImageType  > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(colorStrokesFilename);
  writer->SetInput(colorStrokeImage);
  writer->Update();

}


void Form::btnPropagate_clicked()
{
  // Give the data to the structure propagation algorithm
  this->StructurePropagationFilter->SetMask(this->Mask);
  this->StructurePropagationFilter->SetPatchRadius(3);
  this->StructurePropagationFilter->SetPropagationLine(this->ColorPropagationLine);

  this->StructurePropagationFilter->PropagateStructure();
  //this->StructurePropagationFilter->PropagateStructure<ColorImageType>();

}

void Form::actionFlip_Image_triggered()
{
  double up[3];
  this->LeftRenderer->GetActiveCamera()->GetViewUp(up);
  up[1] *= -1;
  this->LeftRenderer->GetActiveCamera()->SetViewUp(up);
  this->RightRenderer->GetActiveCamera()->SetViewUp(up);

  double pos[3];
  this->LeftRenderer->GetActiveCamera()->GetPosition(pos);
  pos[2] *= -1;
  this->LeftRenderer->GetActiveCamera()->SetPosition(pos);
  this->RightRenderer->GetActiveCamera()->SetPosition(pos);

  this->Refresh();

  /*
  this->CameraUp[1] *= -1;
  std::cout << this->CameraUp[0] << " " << this->CameraUp[1] << " " << this->CameraUp[2] << std::endl;

  double pos[3];
  this->LeftRenderer->GetActiveCamera()->GetPosition(pos);
  std::cout << "position: " << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
  pos[2] *= -1;
  this->LeftRenderer->GetActiveCamera()->SetPosition(pos);
  this->RightRenderer->GetActiveCamera()->SetPosition(pos);

  this->LeftRenderer->GetActiveCamera()->SetViewUp(this->CameraUp);
  this->RightRenderer->GetActiveCamera()->SetViewUp(this->CameraUp);
  this->Refresh();
  */
}

void Form::actionSave_Result_triggered()
{
  // Ask the user for a filename to save the segment mask image to

  QString fileName = QFileDialog::getSaveFileName(this,
    tr("Save Segment Mask Image"), "/home/doriad", tr("Image Files (*.png *.bmp)"));
/*
  // Convert the image from a 1D vector image to an unsigned char image
  typedef itk::CastImageFilter< GrayscaleImageType, itk::Image<itk::CovariantVector<unsigned char, 1>, 2 > > CastFilterType;
  CastFilterType::Pointer castFilter = CastFilterType::New();
  castFilter->SetInput(this->GraphCut->GetSegmentMask());

  typedef itk::NthElementImageAdaptor< itk::Image<itk:: CovariantVector<unsigned char, 1>, 2 >,
    unsigned char> ImageAdaptorType;

  ImageAdaptorType::Pointer adaptor = ImageAdaptorType::New();
  adaptor->SelectNthElement(0);
  adaptor->SetImage(castFilter->GetOutput());
*/

/*
  // Write the file
  //typedef  itk::ImageFileWriter< ImageAdaptorType > WriterType;
  typedef  itk::ImageFileWriter< MaskImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(fileName.toStdString());
  //writer->SetInput(adaptor);
  writer->SetInput(this->GraphCut->GetSegmentMask());
  writer->Update();
*/
}

void Form::actionOpen_Grayscale_Image_triggered()
{
  // Instantiate everything to use grayscale images and grayscale histograms
  OpenFile<GrayscaleImageType>();
}

void Form::actionOpen_Color_Image_triggered()
{
  // Instantiate everything to use color images and color histograms
  OpenFile<ColorImageType>();
}

template <typename TImageType>
void Form::OpenFile()
{
  // Get a filename to open
  QString filename = QFileDialog::getOpenFileName(this,
     //tr("Open Image"), "/media/portable/Projects/src/InteractiveImageGraphCutSegmentation/data", tr("Image Files (*.png *.bmp *.mhd)"));
      tr("Open Image"), ".", tr("Image Files (*.png *.bmp *.mhd)"));

  if(filename.isEmpty())
    {
    return;
    }

  // Clear the scribbles
  //this->ScribbleInteractorStyle->ClearStrokes();

  // Read file
  typename itk::ImageFileReader<TImageType>::Pointer reader = itk::ImageFileReader<TImageType>::New();

  reader->SetFileName(filename.toStdString());
  reader->Update();

  this->ImageRegion = reader->GetOutput()->GetLargestPossibleRegion();

  /*
  // Delete the old object if one exists
  if(this->GraphCut)
    {
    delete this->GraphCut;
    }

  // Instantiate the ImageGraphCut object with the correct type
  this->GraphCut = new ImageGraphCut<TImageType>(reader->GetOutput());
  */
  // Convert the ITK image to a VTK image and display it
  vtkSmartPointer<vtkImageData> VTKImage =
    vtkSmartPointer<vtkImageData>::New();
  Helpers::ITKImageToVTKImage<TImageType>(reader->GetOutput(), VTKImage);

  this->LeftRenderer->RemoveAllViewProps();

  this->OriginalImageActor->SetInput(VTKImage);
  this->ScribbleInteractorStyle->InitializeTracer(this->OriginalImageActor);

  this->LeftRenderer->ResetCamera();
  this->Refresh();

  this->StructurePropagationFilter = new StructurePropagation<TImageType>(reader->GetOutput());

}

void Form::Refresh()
{
  this->LeftRenderer->Render();
  this->RightRenderer->Render();
  this->qvtkWidgetRight->GetRenderWindow()->Render();
  this->qvtkWidgetLeft->GetRenderWindow()->Render();
  this->qvtkWidgetRight->GetInteractor()->Render();
  this->qvtkWidgetLeft->GetInteractor()->Render();
}

void Form::StrokeUpdated(vtkPolyData* path, bool closed)
{
  std::cout << "Form::StrokeUpdated" << std::endl;
  if(this->radDrawHole->isChecked())
    {
    UpdateMaskFromStroke(path, closed);
    }
  else if(this->radDrawPropagationLine->isChecked())
    {
    UpdateColorPropagationLineFromStroke(path);
    }
}


void Form::UpdateColorPropagationLineFromStroke(vtkPolyData* polyDataPath)
{
  // For now, we only store one path at a time. In the future, we'd want use a std::vector<std::vector<itk::Index<2> > > paths
  // and do this->Paths.push_back(path);

  //this->ColorPropagationPathActor->SetInput(polyDataPath);
  this->ColorPropagationPathPolyData->ShallowCopy(polyDataPath);
  std::cout << this->ColorPropagationPathPolyData->GetNumberOfPoints() << " points." << std::endl;

  std::vector<itk::Index<2> > path = Helpers::PolyDataToPixelList(polyDataPath);
  this->ColorPropagationLine = path;

  this->LeftRenderer->AddActor(this->ColorPropagationPathActor);
  this->Refresh();
}

void Form::UpdateMaskFromStroke(vtkPolyData* path, bool closed)
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