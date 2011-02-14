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

template <typename TImage>
void ComputationThread<TImage>::run()
{
  // When the thread is started, emit the signal to start the marquee progress bar
  //std::cout << "emitting StartProgressSignal..." << std::endl;
  Q_EMIT StartProgressSignal();

  this->Propagation->PropagateStructure();

  // When the function is finished, end the thread
  exit();
}

template <typename TImage>
void ComputationThread<TImage>::exit()
{
  // When the thread is stopped, emit the signal to stop the marquee progress bar
  std::cout << "emitting StopProgressSignal..." << std::endl;
  Q_EMIT StopProgressSignal();
}
