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

#include "fflComputeN2.cuh"

#define RX (xSh[threadIdx.x] - xShCh[i])
#define RY (ySh[threadIdx.x] - yShCh[i])
#define RZ (zSh[threadIdx.x] - zShCh[i])

#define IND (blockDim.x * blockIdx.x + threadIdx.x)

#define SQ(x) ((x)*(x))

__constant__ ConstParams devCPar;
__constant__ ChangableParams devChPar;

const float myu0 = 4e-7 * M_PI;
const float g = 9.81;

__device__ int hCounter = 0;


__shared__ float xSh[SHARED_ARRAY];
__shared__ float ySh[SHARED_ARRAY];
__shared__ float zSh[SHARED_ARRAY];

__shared__ float x1Sh[SHARED_ARRAY];
__shared__ float y1Sh[SHARED_ARRAY];
__shared__ float z1Sh[SHARED_ARRAY];

__shared__ float uxSh[SHARED_ARRAY];
__shared__ float uySh[SHARED_ARRAY];
__shared__ float uzSh[SHARED_ARRAY];

__shared__ float xShCh[SHARED_ARRAY];
__shared__ float yShCh[SHARED_ARRAY];
__shared__ float zShCh[SHARED_ARRAY];

__shared__ float x1ShCh[SHARED_ARRAY];
__shared__ float y1ShCh[SHARED_ARRAY];
__shared__ float z1ShCh[SHARED_ARRAY];

__shared__ float uxShCh[SHARED_ARRAY];
__shared__ float uyShCh[SHARED_ARRAY];
__shared__ float uzShCh[SHARED_ARRAY];


__global__ void setupKernel(curandState *state, unsigned long seed) {
	register int ind;
	ind = IND;
	curand_init (seed, ind, 0, &state[ind]);
}



__device__ fh getFH(PartParams mtx) {
	fh ret;
/*	if (devCPar.gravitation == true) {
		ret.f.z = - devCPar.r * g * (devCPar.roParticles - devCPar.roEnvironment)
				/ (4.0f * M_PI * devCPar.myu * devCPar.myu * myu0 / 3.0f);
	}*/

	register float b, b2, d1, d2, dd5, dd7;
	register float bb2d7, umd5;
	register float pow3s2d2, flj;
	register float rx, ry, rz;
	register int ind;
	ind = IND;

	for (unsigned int k = 0; k < devCPar.nPart; k += blockDim.x ) {
		xShCh[threadIdx.x] = mtx.x[threadIdx.x + k];
		yShCh[threadIdx.x] = mtx.y[threadIdx.x + k];
		zShCh[threadIdx.x] = mtx.z[threadIdx.x + k];

		uxShCh[threadIdx.x] = mtx.ux[threadIdx.x + k];
		uyShCh[threadIdx.x] = mtx.uy[threadIdx.x + k];
		uzShCh[threadIdx.x] = mtx.uz[threadIdx.x + k];
		__syncthreads();

		for (unsigned int i = 0; i < blockDim.x; i++) {
			if (k + i != ind) {
				rx = xSh[threadIdx.x] - xShCh[i];
				ry = ySh[threadIdx.x] - yShCh[i];
				rz = zSh[threadIdx.x] - zShCh[i];
				b = uxShCh[i] * rx + uyShCh[i] * ry + uzShCh[i] * rz;
				b2 = uxSh[threadIdx.x] * rx + uySh[threadIdx.x] * ry + uzSh[threadIdx.x] * rz;

				d2 = rx * rx + ry * ry + rz * rz;
				d1 = __powf(d2, 0.5f);

				dd5 = __fdividef(1.0f, d2 * d2 * d1);
				dd7 = __fdividef(dd5, d2);

				bb2d7 = 15.0f * b * b2 * dd7;

				umd5 = - 3.0f * dd5
						* (uxShCh[i] * uxSh[threadIdx.x] + uyShCh[i] * uySh[threadIdx.x] + uzShCh[i] * uzSh[threadIdx.x]);

				ret.h.x += (b * 3.0f * rx - d2 * uxShCh[i]) * dd5;
				ret.h.y += (b * 3.0f * ry - d2 * uyShCh[i]) * dd5;
				ret.h.z += (b * 3.0f * rz - d2 * uzShCh[i]) * dd5;

				ret.f.x += -rx * (umd5 + bb2d7)
						+ 3.0f * (b * uxSh[threadIdx.x] + b2 * uxShCh[i]) * dd5;

				ret.f.y += -ry * (umd5  +  bb2d7)
						+ 3.0f * (b * uySh[threadIdx.x] + b2 * uyShCh[i]) * dd5;

				ret.f.z += -rz * (umd5  +  bb2d7)
						+ 3.0f * (b * uzSh[threadIdx.x] + b2 * uzShCh[i]) * dd5;

				if (d1 < devCPar.ljCutOffR * devCPar.sigma) {
					/* Lennard Jonnes force is equal:
					 * flj = 24.0f * E * (2 * __powf(devCPar.sigma2/D2, 6) - __powf(devCPar.sigma2/D2, 3)) / d1;
					 * where devCPar.sigma2 is devCPar.sigma*devCPar.sigma;
					 * But each component is equal to (for x) fx = RX  * flj / D1;
					 * the last D1 in were included to force flj to decrease calculations;
					 * */
					pow3s2d2 = __powf(__fdividef(SQ(devCPar.sigma), d2), 3);

					flj = __fdividef(24.0f * devCPar.eps * (pow3s2d2  * (2.0f * pow3s2d2 - 1.0f)), d2);
					flj -= devCPar.ljCutOffForce;
					ret.f.x +=  rx * flj;
					ret.f.y +=  ry * flj;
					ret.f.z +=  rz * flj;
				}
			}
		}
		__syncthreads();
	}

	if (ind == 0 && hCounter > 99) {
		printf("hx = %f, hy = %f, hz = %f\n", ret.f.x, ret.f.y, ret.f.z);
		hCounter = 0;
	}

	if (ind == 0) {
		hCounter++;
	}

	ret.h.x += devChPar.hExtX;
	ret.h.y += devChPar.hExtY;
	ret.h.z += devChPar.hExtZ;

	if(devCPar.cf == ConstParams::BARREL) {
		//making center (0;0)
		xSh[threadIdx.x] = xSh[threadIdx.x] - devCPar.barrelR;
		ySh[threadIdx.x] = ySh[threadIdx.x] - devCPar.barrelR;
		float xr, yr;
		//calculating nearest point of the circle to current particle
		xr = __powf(SQ(devCPar.barrelR) /
				(1.0f + SQ(ySh[threadIdx.x]/(xSh[threadIdx.x] + 1E-8))), 0.5f);
		yr = xr * ySh[threadIdx.x] / (xSh[threadIdx.x] + 1E-8);

		//there are 2 roots, we should check + and -
		if (SQ(xSh[threadIdx.x] - xr) + SQ(ySh[threadIdx.x] - yr)
				> SQ(xSh[threadIdx.x] + xr) + SQ(ySh[threadIdx.x] + yr)) {
			xr = -xr;
			yr = -yr;
		}

		rx = xSh[threadIdx.x] - xr;
		ry = ySh[threadIdx.x] - yr;
		d2 = (rx * rx + ry * ry);

		pow3s2d2 = __powf(SQ(devCPar.sigmaWall)/d2, 3);
		flj = 24.0f * devCPar.epsWall * (pow3s2d2  * (2.0f * pow3s2d2 - 1.0f)) / d2;
		ret.f.x +=  rx * flj;
		ret.f.y +=  ry * flj;

		ret.f.z += 24.0f * devCPar.epsWall * (2.0f * __powf(devCPar.sigmaWall/zSh[threadIdx.x], 12) - __powf(devCPar.sigmaWall/zSh[threadIdx.x], 6)) / zSh[threadIdx.x];
		ret.f.z += -24.0f * devCPar.epsWall * (2.0f * __powf(devCPar.sigmaWall/(devCPar.lz - zSh[threadIdx.x]), 12) - __powf(devCPar.sigmaWall/(devCPar.lz - zSh[threadIdx.x]), 6)) / (devCPar.lz -zSh[threadIdx.x]);
	} else {
		ret.f.x += 24.0f * devCPar.epsWall * (2.0f * __powf(devCPar.sigmaWall/xSh[threadIdx.x], 12) - __powf(devCPar.sigmaWall/xSh[threadIdx.x], 6)) / xSh[threadIdx.x];
		ret.f.y += 24.0f * devCPar.epsWall * (2.0f * __powf(devCPar.sigmaWall/ySh[threadIdx.x], 12) - __powf(devCPar.sigmaWall/ySh[threadIdx.x], 6)) / ySh[threadIdx.x];
		ret.f.z += 24.0f * devCPar.epsWall * (2.0f * __powf(devCPar.sigmaWall/zSh[threadIdx.x], 12) - __powf(devCPar.sigmaWall/zSh[threadIdx.x], 6)) / zSh[threadIdx.x];

		ret.f.x += -24.0f * devCPar.epsWall * (2.0f * __powf(devCPar.sigmaWall/(devCPar.lx - xSh[threadIdx.x]), 12) - __powf(devCPar.sigmaWall/(devCPar.lx - xSh[threadIdx.x]), 6)) / (devCPar.lx -xSh[threadIdx.x]);
		ret.f.y += -24.0f * devCPar.epsWall * (2.0f * __powf(devCPar.sigmaWall/(devCPar.ly - ySh[threadIdx.x]), 12) - __powf(devCPar.sigmaWall/(devCPar.ly - ySh[threadIdx.x]), 6)) / (devCPar.ly -ySh[threadIdx.x]);
		ret.f.z += -24.0f * devCPar.epsWall * (2.0f * __powf(devCPar.sigmaWall/(devCPar.lz - zSh[threadIdx.x]), 12) - __powf(devCPar.sigmaWall/(devCPar.lz - zSh[threadIdx.x]), 6)) / (devCPar.lz -zSh[threadIdx.x]);
	}
	return ret;
}


__global__
__launch_bounds__(SHARED_ARRAY)
void averageUKernel(PartParams mtx, float3* u) {
	float sinTheta, cosTheta;
	float sinPhy, cosPhy;
	register int ind;
	ind = IND;
	__shared__ float ux[SHARED_ARRAY], uy[SHARED_ARRAY], uz[SHARED_ARRAY];

	if(threadIdx.x + blockIdx.x * blockDim.x == 0) {
		u->x = 0.0f;
		u->y = 0.0f;
		u->z = 0.0f;
	}

	__sincosf(mtx.theta[ind], &sinTheta, &cosTheta);
	__sincosf(mtx.phy[ind], &sinPhy, &cosPhy);

	ux[threadIdx.x] = sinTheta * cosPhy;
	uy[threadIdx.x] = sinTheta * sinPhy;
	uz[threadIdx.x] = cosTheta;

	__syncthreads();

	for (int cntr = SHARED_ARRAY / 2; cntr > 0; cntr /= 2 ) {
		if (threadIdx.x < cntr) {
			if (threadIdx.x  + cntr < blockDim.x) {
//				printf("Sum[%d]: %d and %d \n", blockIdx.x, threadIdx.x, threadIdx.x + cntr);
				ux[threadIdx.x] = 0.5 * (ux[threadIdx.x] + ux[threadIdx.x + cntr]);
				uy[threadIdx.x] = 0.5 * (uy[threadIdx.x] + uy[threadIdx.x + cntr]);
				uz[threadIdx.x] = 0.5 * (uz[threadIdx.x] + uz[threadIdx.x + cntr]);
			}
		}
		__syncthreads();
	}

	if (threadIdx.x == 0) {
		ux[0] /= gridDim.x;
		uy[0] /= gridDim.x;
		uz[0] /= gridDim.x;
		atomicAdd(&(u->x), ux[0]);
		atomicAdd(&(u->y), uy[0]);
		atomicAdd(&(u->z), uz[0]);
	}
}

__global__
__launch_bounds__(SHARED_ARRAY)
void calcU(PartParams mtx) {
	float sinTheta, cosTheta;
	float sinPhy, cosPhy;
	register int ind;
	ind = IND;
	__sincosf(mtx.theta[ind], &sinTheta, &cosTheta);
	__sincosf(mtx.phy[ind], &sinPhy, &cosPhy);
	int id = ind;
	mtx.ux[id] = sinTheta * cosPhy;
	mtx.uy[id] = sinTheta * sinPhy;
	mtx.uz[id] = cosTheta;
}

__global__
__launch_bounds__(SHARED_ARRAY)
void applyDeltas(PartParams mtx) {
	register int ind;
	ind = IND;

	mtx.phy[ind] += mtx.deltaPhy[ind];
	mtx.theta[ind] += mtx.deltaTheta[ind];
	mtx.x[ind] += mtx.deltaX[ind];
	mtx.y[ind] += mtx.deltaY[ind];
	mtx.z[ind] += mtx.deltaZ[ind];
}

__global__
__launch_bounds__(SHARED_ARRAY)
void oneStep(PartParams mtx)
{
	register int ind;
	ind = IND;

	xSh[threadIdx.x] = mtx.x[ind];
	ySh[threadIdx.x] = mtx.y[ind];
	zSh[threadIdx.x] = mtx.z[ind];

	x1Sh[threadIdx.x] = mtx.x1[ind];
	y1Sh[threadIdx.x] = mtx.y1[ind];
	z1Sh[threadIdx.x] = mtx.z1[ind];

	uxSh[threadIdx.x] = mtx.ux[ind];
	uySh[threadIdx.x] = mtx.uy[ind];
	uzSh[threadIdx.x] = mtx.uz[ind];

	register fh infl;
	infl = getFH(mtx);

	register float3 n;

	n.x = uySh[threadIdx.x] * infl.h.z - uzSh[threadIdx.x] * infl.h.y;
	n.y = uzSh[threadIdx.x] * infl.h.x - uxSh[threadIdx.x] * infl.h.z;
	n.z = uxSh[threadIdx.x] * infl.h.y - uySh[threadIdx.x] * infl.h.x;

	mtx.deltaPhy[ind] = mtx.phy1[ind] * devChPar.dTimeCurrent;
	mtx.deltaTheta[ind] = mtx.theta1[ind] * devChPar.dTimeCurrent;

	mtx.deltaX[ind] = x1Sh[threadIdx.x] * devChPar.dTimeCurrent;
	mtx.deltaY[ind] = y1Sh[threadIdx.x] * devChPar.dTimeCurrent;
	mtx.deltaZ[ind] = z1Sh[threadIdx.x] * devChPar.dTimeCurrent;

	if (devCPar.thermalBath == true) {
		curandState localState = mtx.rndStates[ind];

		mtx.phy1[ind] += 2.5f * ((n.z - mtx.phy1[ind] * devCPar.nyu) * devChPar.dTimeCurrent
				+ curand_normal(&localState) /* * 3.0f * devCPar.nyu*/ * devChPar.sqrtdTime * devCPar.qr);
		mtx.theta1[ind] += 2.5f * ((-n.x * __sinf(mtx.phy[ind]) + n.y * __cosf(mtx.phy[ind])
					- mtx.theta1[ind] * devCPar.nyu) * devChPar.dTimeCurrent
					+ curand_normal(&localState) /* * 3.0f * devCPar.nyu*/ * devChPar.sqrtdTime * devCPar.qr);

		mtx.x1[ind] += (infl.f.x - x1Sh[threadIdx.x] * devCPar.eta) * devChPar.dTimeCurrent
				+ curand_normal(&localState) /* * 2.0f * devCPar.eta*/ * devChPar.sqrtdTime * devCPar.qt;
		mtx.y1[ind] += (infl.f.y - y1Sh[threadIdx.x] * devCPar.eta) * devChPar.dTimeCurrent
				+ curand_normal(&localState) /* * 2.0f * devCPar.eta*/ * devChPar.sqrtdTime * devCPar.qt;
		mtx.z1[ind] += (infl.f.z - z1Sh[threadIdx.x] * devCPar.eta) * devChPar.dTimeCurrent
				+ curand_normal(&localState) /* * 2.0f * devCPar.eta*/ * devChPar.sqrtdTime * devCPar.qt;
		mtx.rndStates[ind] = localState;
	} else {
		mtx.phy1[ind] += 2.5f * (n.z - mtx.phy1[ind] * devCPar.nyu) * devChPar.dTimeCurrent;

		mtx.theta1[ind] += 2.5f * (-n.x * __sinf(mtx.phy[ind]) + n.y * __cosf(mtx.phy[ind])
					- mtx.theta1[ind] * devCPar.nyu) * devChPar.dTimeCurrent;

		mtx.x1[ind] += (infl.f.x - x1Sh[threadIdx.x] * devCPar.eta) * devChPar.dTimeCurrent;
		mtx.y1[ind] += (infl.f.y - y1Sh[threadIdx.x] * devCPar.eta) * devChPar.dTimeCurrent;
		mtx.z1[ind] += (infl.f.z - z1Sh[threadIdx.x] * devCPar.eta) * devChPar.dTimeCurrent;
	}
}


__global__
__launch_bounds__(SHARED_ARRAY)
void uploadSharedMemory(PartParams mtx) {
	register int ind;
	ind = IND;

	xSh[threadIdx.x] = mtx.x[ind];
	ySh[threadIdx.x] = mtx.y[ind];
	zSh[threadIdx.x] = mtx.z[ind];

	x1Sh[threadIdx.x] = mtx.x1[ind];
	y1Sh[threadIdx.x] = mtx.y1[ind];
	z1Sh[threadIdx.x] = mtx.z1[ind];

	uxSh[threadIdx.x] = mtx.ux[ind];
	uySh[threadIdx.x] = mtx.uy[ind];
	uzSh[threadIdx.x] = mtx.uz[ind];
}

__global__
__launch_bounds__(SHARED_ARRAY)
void getCurrentPhysKernel(PartParams mtx, float* ringPhy){
	float x, y, phy;
	int ind;
	ind = IND;
	x = mtx.x[ind] - devCPar.barrelR;
	y = mtx.y[ind] - devCPar.barrelR;
	phy = atan2f(y, x);
	ringPhy[ind] = phy;
}

__global__
__launch_bounds__(SHARED_ARRAY)
void getRingDistrKernel(PartParams mtx, float* oldPhy, float* newPhy, int nRings, float* dphy, int* partPos) {
	float dPhyCached, r;
	int interval;
	int ind = IND;
	__shared__ float dPhyAddayShared[SHARED_ARRAY];	//nRings should not be bigger than SHARED_ARRAY;
	__shared__ int sharedCounter[SHARED_ARRAY];

	if (threadIdx.x < nRings) {
		dPhyAddayShared[threadIdx.x] = 0.0f;
		sharedCounter[threadIdx.x] = 0;
	}
	if (ind < nRings) {
		dphy[ind] = 0.0f;
		partPos[ind] = 0;
	}

	xSh[threadIdx.x] = mtx.x[ind];
	ySh[threadIdx.x] = mtx.y[ind];
	__syncthreads();

	dPhyCached = newPhy[ind] - oldPhy[ind];
	if (dPhyCached > M_PI) dPhyCached -= 2.0f * M_PI;
	if (dPhyCached < -M_PI) dPhyCached += 2.0f * M_PI;

	r = __powf(SQ(xSh[threadIdx.x] - devCPar.barrelR) + SQ(ySh[threadIdx.x] - devCPar.barrelR), 0.5f);
	interval = (int) floorf((float)nRings * r / devCPar.barrelR);	//ly - is R of barrel
	atomicAdd(&(dPhyAddayShared[interval]), dPhyCached);
	atomicAdd(&(sharedCounter[interval]), 1);
	__syncthreads();

	if (threadIdx.x == 0) {
		for (int i = 0; i < nRings; i++) {
			if (sharedCounter[i] > 0) atomicAdd(&(dphy[i]), dPhyAddayShared[i] / (sharedCounter[i] * gridDim.x));
			atomicAdd(&(partPos[i]), sharedCounter[i]);
		}
	}

}

void fillGloabalChangable(ChangableParams* chPar) {
	cudaMemcpyToSymbol(devChPar, chPar, sizeof(ChangableParams), 0, cudaMemcpyHostToDevice);
}

void fillGloabalConstant(ConstParams* cPar) {
	cudaMemcpyToSymbol(devCPar, cPar, sizeof(ConstParams), 0, cudaMemcpyHostToDevice);
}

void runCalcU(int blocks, int threads, PartParams devMatrixes) {
	calcU<<<blocks, threads>>>(devMatrixes);
}
void runOneStep(int blocks, int threads, PartParams devMatrixes) {
	oneStep<<<blocks, threads>>>(devMatrixes);
}
void runApplyDeltas(int blocks, int threads, PartParams devMatrixes) {
	applyDeltas<<<blocks, threads>>>(devMatrixes);
}

void runSetupKernel(int blocks, int threads, PartParams devMatrixes) {
	setupKernel<<<blocks, threads>>>(devMatrixes.rndStates, time(NULL));
}


float3 averageU(int blocks, int threads, PartParams devMatrixes) {
	float3 u;
	float3* devU;
	cudaMalloc((void**)&devU, sizeof(float3));
	averageUKernel<<<blocks, threads>>>(devMatrixes, devU);
	cudaMemcpy(&u, devU, sizeof(float3), cudaMemcpyDeviceToHost);
	cudaFree(devU);
	return u;
}

void getCurrentPhys(int blocks, int threads, PartParams devMatrixes, float** devOldPhy) {
	getCurrentPhysKernel<<<blocks, threads>>> (devMatrixes, *devOldPhy);
}

void getRingStat(int blocks, int threads, PartParams devMatrixes, int nRings, int nPart, float** devOldPhy, float* dphy, int* pd) {
	float* devNewPhy;
	float* devDPhy;
	int* devPD;
	float* tmp;

	if (cudaSuccess != cudaMalloc((void**)&devNewPhy, sizeof(float) * nPart))
		throw DeviceMemoryAllocationException("Error allocation of devNewPhy");
	if (cudaSuccess != cudaMalloc((void**)&devDPhy, sizeof(float) * nRings))
		throw DeviceMemoryAllocationException("Error allocation of devDPhy");
	if (cudaSuccess != cudaMalloc((void**)&devPD, sizeof(int) * nRings))
		throw DeviceMemoryAllocationException("Error allocation of devPD");

	getCurrentPhys(blocks, threads, devMatrixes, &devNewPhy);
	getRingDistrKernel<<<blocks, threads>>>(devMatrixes, *devOldPhy, devNewPhy, nRings, devDPhy, devPD);

	cudaMemcpy(dphy, devDPhy, sizeof(float) * nRings, cudaMemcpyDeviceToHost);
	cudaMemcpy(pd, devPD, sizeof(int) * nRings, cudaMemcpyDeviceToHost);

	tmp = *devOldPhy;
	*devOldPhy = devNewPhy;


	if (cudaSuccess != cudaFree(tmp))
		throw DeviceMemoryException(std::string("Error free of tmp") +
				std::string(cudaGetErrorString(cudaGetLastError())));
	if (cudaSuccess != cudaFree(devDPhy))
		throw DeviceMemoryException(std::string("Error free of devDPhy") +
				std::string(cudaGetErrorString(cudaGetLastError())));
	if (cudaSuccess != cudaFree(devPD))
		throw DeviceMemoryException(std::string("Error free of devPD") +
				std::string(cudaGetErrorString(cudaGetLastError())));
}
