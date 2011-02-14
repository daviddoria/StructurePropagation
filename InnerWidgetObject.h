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
 * This class exists because Qt does not allow templated classes to declare signals or slots.
 * We define pure virtual slots which are defined in InnerWidget.
 */

#ifndef INNERWIDGETOBJECT_H
#define INNERWIDGETOBJECT_H

// Qt
#include "ui_InnerWidget.h"
#include <QWidget>

class InnerWidgetObject : public QWidget, public Ui::InnerWidget
{
Q_OBJECT
public:
  InnerWidgetObject(QWidget *parent = 0);
  ~InnerWidgetObject(){}

public Q_SLOTS:

  // Buttons
  virtual void btnClearStrokes_clicked() = 0;
  virtual void btnSaveStrokes_clicked() = 0;
  virtual void btnPropagate_clicked() = 0;
  virtual void btnLoadMask_clicked() = 0;
  virtual void btnSaveResult_clicked() = 0;
  virtual void btnScreenshot_clicked() = 0;

  // Other GUI elements
  virtual void sldPatchRadius_valueChanged() = 0;

  // Checkboxes
  virtual void chkFlip_clicked() = 0;
  virtual void chkScale_clicked() = 0;
  virtual void chkShowPaths_clicked() = 0;
  virtual void chkShowResult_clicked() = 0;
  virtual void chkShowOriginal_clicked() = 0;
  virtual void chkShowMask_clicked() = 0;

  // These slots handle running the progress bar while the computations are done in a separate thread.
  virtual void StartProgressSlot() = 0;
  virtual void StopProgressSlot() = 0;

};

#endif