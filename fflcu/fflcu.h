/*
 * fflcu.h
 *
 *  Created on: Jan 12, 2012
 *      Author: alexander
 */

#ifndef FFLCU_H_
#define FFLCU_H_

#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <curand.h>
#include <curand_kernel.h>
#include <time.h>

#include "ConfigFile.h"
#include "fflComputeN2.cuh"
#include "fflComputeBH.cuh"
#include "Exceptions.h"


class Fflcun2 {
public:
	struct ParticlesPosition{
		float *x;
		float *y;
		float *z;
		float *theta;
		float *phy;
		int n;
	};


	Fflcun2(std::string fname);
	Fflcun2();
	virtual ~Fflcun2();
	virtual ParticlesPosition run(int nSteps);
	virtual void setConfig(std::string fname);

	void setHx(float hx) {this->chPar->hExtX = hx;};
	void setHy(float hy) {this->chPar->hExtY = hy;};
	void setHz(float hz) {this->chPar->hExtZ = hz;};
	void setEta(float eta) {this->cPar->eta = eta;};
	void setNyu(float nyu) {this->cPar->nyu = nyu;};
	void setLx(float lx) {this->cPar->lx = lx;};
	void setLy(float ly) {this->cPar->ly = ly;};
	void setLz(float lz) {this->cPar->lz = lz;};
	void setNRings(int nRings) {this->nRings = nRings;};
	virtual void initRotStat();
	void destroyRotStat();
	virtual void getRotStat(float* dPhy, int* pos);
	virtual void setN(int n);
	ParticlesPosition getState() {return this->partPos;}
	void saveConf(std::fstream& out);
	virtual float3 getAverageU();

protected:
	bool checkParticlesLocation(unsigned int maxN);
	virtual void initDevMatrixes();
	virtual void allocMemory();
	virtual void freeMemory();
	virtual void setSettings();

	ConstParams *cPar;
	ChangableParams *chPar;

	PartParams devMatrixes;
	int blocks, threads;
	int attempts;

	ParticlesPosition partPos;
	std::string configFileName;

	float* devOldPhy;
	int nRings;
	int deviceNubmer;
	double carrierViscosity;
};

class Fflcubh: public Fflcun2 {
public:
	Fflcubh(std::string fname);
	/*virtual*/ void setN(int n);
	/*virtual*/ ~Fflcubh();
	/*virtual*/ ParticlesPosition run(int nSteps);
	/*virtual*/ float3 getAverageU();
	/*virtual*/ void initRotStat();
	/*virtual*/ void getRotStat(float* dPhy, int* pos);

private:
	/*virtual*/ void allocMemory();
	/*virtual*/ void freeMemory();
	/*virtual*/ void initDevMatrixes();
	Box boxl;
	BHArrays arrl;
	int nnodes;
	float *massl;
};

#endif /* FFLCU_H_ */
