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
