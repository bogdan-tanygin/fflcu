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
 * DataTypes.h
 *
 *  Created on: 24 окт. 2011
 *      Author: alexander
 */

#ifndef DATATYPES_H_
#define DATATYPES_H_

#include <curand_kernel.h>

struct ConstParams{
	enum ContainerForm {BRICK, BARREL };
	enum FieldType {PERMANENT, ROTATING};
	int nPart;	//number of particles

	float dTimeInitial;
	float timeMax;

	float lx;
	float ly;
	float lz;

	float nyu;
	float eta;	//0.1

	float r;
	float myu;
	float roParticles;
	float roEnvironment;

	bool gravitation;

	bool thermalBath;

	ContainerForm cf;
	float barrelR;

	float sigma;
	float sigmaWall;
	float eps;	//Lennard Jones particles energy
	float epsWall;	//Lennard jones wall energy
	FieldType ft;	//field type
	float rff;	//rotating field frequency

	float hExtXInit;
	float hExtYInit;
	float hExtZInit;

	float ljCutOffR;
	float ljCutOffForce;
	float ljCutOffForceWall;
	//float factorTime;

	float temper;	//temperature in K

	float qr;	//noise parameter
	float qt;
};

struct ChangableParams{
	float time;
	float dTimeCurrent;
	float sqrtdTime;
	float hExtX;
	float hExtY;
	float hExtZ;
};


struct PartParams{
	float* x;
	float* y;
	float* z;
	float* phy;
	float* theta;
	float* ux;
	float* uy;
	float* uz;

	float* phy1; //=d(phy)/(dt)
	float* theta1; //=d(theta)/(dt)
	float* x1;
	float* y1;
	float* z1;

	float* deltaX;
	float* deltaY;
	float* deltaZ;
	float* deltaPhy;
	float* deltaTheta;

	curandState* rndStates;
};

#endif /* DATATYPES_H_ */

