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
  // Menu items
  //virtual void actionFlip_Image_triggered() = 0;
  //virtual void actionSave_Result_triggered() = 0;

  // Buttons
  virtual void btnClearStrokes_clicked() = 0;
  virtual void btnSaveStrokes_clicked() = 0;
  virtual void btnPropagate_clicked() = 0;
  virtual void btnLoadMask_clicked() = 0;

  // Radio buttons
  //void radDrawPropagationLine_clicked();
  //void radDrawHole_clicked();

  // These slots handle running the progress bar while the computations are done in a separate thread.
  //virtual void StartProgressSlot() = 0;
  //virtual void StopProgressSlot() = 0;

};

#endif