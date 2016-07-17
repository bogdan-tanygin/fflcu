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
 * Visualization.cpp
 *
 *  Created on: 12 сент. 2011
 *      Author: alexander
 */

#include "ConsoleVisualization.h"

ConsoleVisualization::ConsoleVisualization() {
	return;
}

/*virtual*/ void ConsoleVisualization::showThis(float* x, float* y, float* z, float* theta, float* phy, int n){
	cout<<"{";

	for(int i = 0; i < n; i++){
		cout<< setiosflags(ios::fixed) << setprecision(1) << x[i]<<";"<<y[i]<<";"<<z[i]<<";"<<theta[i]<<";"<<phy[i]<<";";
	}
	cout<<"}"<<endl;
}

/*virtual*/ void ConsoleVisualization::setSettings(int xScr, int yScr, float xSize, float ySize, float zSize){
	return;	//empty function
}

/*virtual*/ void ConsoleVisualization::showCurrentState(float t, int it, float var, float kx, float ky, float kz){
	cout << setiosflags(ios::fixed) << setprecision(2)
			<< "t:\t" << t <<
			"\nIteration:\t" << it <<
			"\nVariable:\t" << var <<
			"\nX-Susceptibility:\t" << kx <<
			"\nY-Susceptibility:\t" << ky <<
			"\nZ-Susceptibility:\t" << kz <<
			"\n=============================" << endl;
}

ConsoleVisualization::~ConsoleVisualization() {

}
