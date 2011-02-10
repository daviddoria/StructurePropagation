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

#include "MainWindow.h"
#include "InnerWidget.h"

MainWindow::MainWindow(QWidget *parent)
{
  // Setup the GUI and connect all of the signals and slots
  setupUi(this);
  connect( this->actionOpen_Color_Image, SIGNAL( triggered() ), this, SLOT(actionOpen_Color_Image_triggered()) );
  connect( this->actionOpen_Grayscale_Image, SIGNAL( triggered() ), this, SLOT(actionOpen_Grayscale_Image_triggered()) );

}

void MainWindow::actionOpen_Grayscale_Image_triggered()
{
  InnerWidget<GrayscaleImageType>* innerWidget = new InnerWidget<GrayscaleImageType>(this);
  this->setCentralWidget(innerWidget);
  innerWidget->OpenFile();
}

void MainWindow::actionOpen_Color_Image_triggered()
{
  InnerWidget<ColorImageType>* innerWidget = new InnerWidget<ColorImageType>(this);
  this->setCentralWidget(innerWidget);
  innerWidget->OpenFile();
}
