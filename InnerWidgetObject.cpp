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

#include "InnerWidgetObject.h"

InnerWidgetObject::InnerWidgetObject(QWidget *parent)
{
  // Setup the GUI and connect all of the signals and slots
  setupUi(this);

  // Buttons
  connect( this->btnPropagate, SIGNAL( clicked() ), this, SLOT(btnPropagate_clicked()));
  connect( this->btnLoadMask, SIGNAL( clicked() ), this, SLOT(btnLoadMask_clicked()));
  connect( this->btnClearStrokes, SIGNAL( clicked() ), this, SLOT(btnClearStrokes_clicked()));
  connect( this->btnSaveStrokes, SIGNAL( clicked() ), this, SLOT(btnSaveStrokes_clicked()));

  // Menu
  //connect( this->actionFlip_Image, SIGNAL( triggered() ), this, SLOT(actionFlip_Image_triggered()));
  //connect( this->actionSave_Result, SIGNAL( triggered() ), this, SLOT(actionSave_Result_triggered()));

  // Thread
  //connect(&ProgressThread, SIGNAL(StartProgressSignal()), this, SLOT(StartProgressSlot()), Qt::QueuedConnection);
  //connect(&ProgressThread, SIGNAL(StopProgressSignal()), this, SLOT(StopProgressSlot()), Qt::QueuedConnection);

}