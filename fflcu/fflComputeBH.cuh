/*
 * fflComputeBH.cuh
 *
 *  Created on: Jan 31, 2012
 *      Author: alexander
 */

#ifndef FFLCOMPUTEBH_CUH_
#define FFLCOMPUTEBH_CUH_

/**
 *
 */

#include "DataTypes.cuh"
#include "Exceptions.h"
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <cstdio>
#include <iostream>
#include <time.h>

#define SHARED_ARRAY_BH 512

struct Box{
	float *minx;
	float *maxx;
	float *miny;
	float *maxy;
	float *minz;
	float *maxz;
};

struct BHArrays{
	int *err;
	int *sort;
	int *child;
	int *count;
	int *start;
};

// thread count
#define THREADS0 512	/* must be a power of 2 */
#define THREADS1 512	/* must be a power of 2 */
#define THREADS2 1024
#define THREADS3 1024
#define THREADS4 256
#define THREADS5 128	//256
#define THREADS6 512

// block count = factor * #SMs
#define FACTOR0 1
#define FACTOR1 3
#define FACTOR2 1
#define FACTOR3 1	/* must all be resident at the same time */
#define FACTOR4 1	/* must all be resident at the same time */
#define FACTOR5 5
#define FACTOR6 3

#define WARPSIZE 32
#define MAXDEPTH 32


void fillConstantPointers(int nbodies, int nnodes, float* mass, BHArrays arrl, Box boxl, PartParams mtxl);

void init(int blocks, curandState* devRndStates);

void calcU(int blocks);

void buildBox(int blocks);

void buildTree(int blocks);

void summarize(int blocks);

void sort(int blocks);

void force(int blocks);

void integrate(int blocks);

void fillGloabalChangableBH(ChangableParams* chPar);

void fillGloabalConstantBH(ConstParams* cPar);

float3 averageU(int blocks);

void getCurrentPhysBH(int blocks, float* devOldPhy);

void getRingStatBH(int blocks, int nRings, int nPart, float** devOldPhy, float* dphy, int* pd);

#endif /* FFLCOMPUTEBH_CUH_ */
