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

