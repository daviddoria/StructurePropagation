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

/* The purpose of this class is simply to run the graph cut in a
 * separate thread so that the progress bar can run during the segmentation
*/

#ifndef COMPUTATIONTHREAD_H
#define COMPUTATIONTHREAD_H

#include <QThread>

#include "StructurePropagation.h"

class ComputationThreadObject : public QThread
{
Q_OBJECT
public:

//signals:
Q_SIGNALS:
  // This signal is emitted to start the progress bar
  void StartProgressSignal();

  // This signal is emitted to stop the progress bar
  void StopProgressSignal();
};

template <typename TImage>
class ComputationThread : public ComputationThreadObject
{
public:

  // This function is called when the thread is started
  void run();

  // This function is called when the thread is stopped
  void exit();

  // We need a pointer to this object so we can perform the computations in this thread
  StructurePropagation<TImage>* Propagation;
};

#include "ComputationThread.txx"

#endif