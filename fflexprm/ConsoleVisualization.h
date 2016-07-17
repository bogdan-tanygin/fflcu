/*
  Copyright (C) 2012 Alexander (Polyakov) Peletskyi

  This file is part of FFLCU.

  FFLCU is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  FFLCU is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
 * Visualization.h
 *
 *  Created on: 12 сент. 2011
 *      Author: alexander
 */

#ifndef CONSOLEVISUALIZATION_H_
#define CONSOLEVISUALIZATION_H_

#include <iostream>
#include <iomanip>
#include "IVis.h"

using namespace std;

class ConsoleVisualization: public IVis {
public:
	ConsoleVisualization();
	/*virtual*/ void showThis(float* x, float* y, float* z, float* theta, float* phy, int n);
	/*virtual*/ void setSettings(int xScr, int yScr, float xSize, float ySize, float zSize);
	/*virtual*/ void showCurrentState(float t, int it, float var, float kx, float ky, float kz);
	~ConsoleVisualization();
private:
};

#endif /* CONSOLEVISUALIZATION_H_ */
