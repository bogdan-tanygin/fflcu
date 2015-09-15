/*
 * IVis.h
 *
 *  Created on: Jan 16, 2012
 *      Author: alexander
 */

#ifndef IVIS_H_
#define IVIS_H_

class IVis {
public:
//	virtual IVis() = 0;
	virtual ~IVis(){};
	virtual void showThis(float* x, float* y, float* z, float* theta, float* phy, int n) = 0;
	virtual void showCurrentState(float t, int it, float var, float kx, float ky, float kz) = 0;
	virtual void setSettings(int xScr, int yScr, float xSize, float ySize, float zSize) = 0;
};

#endif /* IVIS_H_ */
