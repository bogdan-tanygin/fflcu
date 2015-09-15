/*
* Copyright 2008, Karen Hains, UWA . All rights reserved.
*
* NOTICE TO USER:
*
* This source code is subject to NVIDIA ownership rights under U.S. and
* international Copyright laws. Users and possessors of this source code
* are hereby granted a nonexclusive, royalty-free license to use this code
* in individual and commercial software.
*
* WE MAKES NO REPRESENTATION ABOUT THE SUITABILITY OF THIS SOURCE
* CODE FOR ANY PURPOSE. IT IS PROVIDED "AS IS" WITHOUT EXPRESS OR
* IMPLIED WARRANTY OF ANY KthreadIdx.x.
*/

#ifndef _FFLCOMPUTE_CUH_
#define _FFLCOMPUTE_CUH_

#include "DataTypes.cuh"
#include "Exceptions.h"
#include <cuda.h>
#include <cstdio>
#include <curand.h>
#include <curand_kernel.h>

#define SHARED_ARRAY 512

struct fh{
	__device__ fh(){
		f.x = 0;
		f.y = 0;
		f.z = 0;
		h.x = 0;
		h.y = 0;
		h.z = 0;
	}
	float3 f;
	float3 h;
};



void fillGloabalChangable(ChangableParams* chPar);
void fillGloabalConstant(ConstParams* cPar);

void runCalcU(int blocks, int threads, PartParams devMatrixes);
void runOneStep(int blocks, int threads, PartParams devMatrixes);
void runApplyDeltas(int blocks, int threads, PartParams devMatrixes);
void runSetupKernel(int blocks, int threads, PartParams devMatrixes);
float3 averageU(int blocks, int threads, PartParams devMatrixes);
void getCurrentPhys(int blocks, int threads, PartParams devMatrixes, float** devOldPhy);
void getRingStat(int blocks, int threads, PartParams devMatrixes, int nRings, int nPart, float** devOldPhy, float* dphy, int* pd);


#endif // #ifndef _FFLCOMPUTE_CUH_
