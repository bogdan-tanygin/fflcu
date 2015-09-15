/*
 * Fflcubh.cpp
 *
 *  Created on: Jan 31, 2012
 *      Author: alexander
 */

#include "fflcu.h"
#include <iostream>
using namespace std;
#define NODES_MULTY 8
//#define DEBUG

Fflcubh::Fflcubh(std::string fname):Fflcun2() {
	int deviceCount;
	cudaDeviceProp deviceProp;

	cudaGetDeviceCount(&deviceCount);
	if (deviceCount == 0) {
		fprintf(stderr, "There is no device supporting CUDA\n");
		exit(-1);
	}

	cudaGetDeviceProperties(&deviceProp, 0);
	if ((deviceProp.major == 9999) && (deviceProp.minor == 9999)) {
		fprintf(stderr, "There is no CUDA capable device\n");
		exit(-1);
	}
	if (deviceProp.major < 2) {
		fprintf(stderr, "Need at least compute capability 2.0\n");
		exit(-1);
	}
	if (deviceProp.warpSize != WARPSIZE) {
		fprintf(stderr, "Warp size must be %d\n", deviceProp.warpSize);
		exit(-1);
	}

	cPar = new ConstParams;
	chPar = new ChangableParams;
	configFileName = fname;
	this->Fflcun2::setSettings();

	blocks = deviceProp.multiProcessorCount;
	fprintf(stderr, "blocks = %d\n", blocks);

	if ((WARPSIZE <= 0) || ((WARPSIZE & (WARPSIZE - 1)) != 0)) {
		fprintf(stderr, "Warp size must be greater than zero and a power of two\n");
		exit(-1);
	}
	if (MAXDEPTH > WARPSIZE) {
		fprintf(stderr, "MAXDEPTH must be less than or equal to WARPSIZE\n");
		exit(-1);
	}
	if ((THREADS1 <= 0) || ((THREADS1 & (THREADS1 - 1)) != 0)) {
		fprintf(stderr, "THREADS1 must be greater than zero and a power of two\n");
		exit(-1);
	}

	nnodes = cPar->nPart * NODES_MULTY;
	if (nnodes < 1024 * blocks) nnodes = 1024 * blocks;
	while ((nnodes & (WARPSIZE - 1)) != 0) nnodes++;
	nnodes--;
	fprintf(stderr, "nnodes = %d, nbodies = %d\n", nnodes, cPar->nPart);

	srand(time(NULL));

	this->allocMemory();
	this->initDevMatrixes();

	fillGloabalChangableBH(chPar);
	fillGloabalConstantBH(cPar);
	fillConstantPointers(cPar->nPart, nnodes, massl, arrl, boxl, devMatrixes);
	init(blocks, devMatrixes.rndStates);
}

Fflcubh::~Fflcubh() {
	cout << "BH destructor" << endl;

	cudaFree(massl);

	cudaFree(arrl.count);
	cudaFree(arrl.sort);
	cudaFree(arrl.start);
	cudaFree(arrl.err);
	cudaFree(arrl.child);

	cudaFree(boxl.maxx);
	cudaFree(boxl.maxy);
	cudaFree(boxl.maxz);
	cudaFree(boxl.minx);
	cudaFree(boxl.miny);
	cudaFree(boxl.minz);
	//freeMemory();
}

void Fflcubh::allocMemory() {
	/* The number of elements in devMatrixes is bigger, than in N2 method.
	 * So, it is necessary to write another method, without using parent one.
	 */
	std::cout << "Memory allocation" << std::endl;
	if (cudaSuccess != cudaMalloc((void **)&arrl.err, sizeof(int)))
		throw DeviceMemoryAllocationException("could not allocate errd\n");	//CudaTest("couldn't allocate errd");
	if (cudaSuccess != cudaMalloc((void **)&arrl.child, sizeof(int) * (nnodes + 1) * 8))
		throw DeviceMemoryAllocationException("could not allocate childd\n");	//CudaTest("couldn't allocate childd");
	if (cudaSuccess != cudaMalloc((void **)&arrl.count, sizeof(int) * (nnodes + 1)))
		throw DeviceMemoryAllocationException("could not allocate countd\n");	//CudaTest("couldn't allocate countd");
	if (cudaSuccess != cudaMalloc((void **)&arrl.start, sizeof(int) * (nnodes + 1)))
		throw DeviceMemoryAllocationException("could not allocate startd\n");	//CudaTest("couldn't allocate startd");
	if (cudaSuccess != cudaMalloc((void **)&arrl.sort, sizeof(int) * (nnodes + 1)))
		throw DeviceMemoryAllocationException("could not allocate sortl\n");	//CudaTest("couldn't allocate sortl");

	if (cudaSuccess != cudaMalloc((void **)&massl, sizeof(float) * (nnodes + 1)))
		throw DeviceMemoryAllocationException("could not allocate massl\n");	//CudaTest("couldn't allocate massl");

	if (cudaSuccess != cudaMalloc((void **)&devMatrixes.x, sizeof(float) * (nnodes + 1)))
		throw DeviceMemoryAllocationException("could not allocate mtxd.x\n");	//CudaTest("couldn't allocate mtxd.x");
	if (cudaSuccess != cudaMalloc((void **)&devMatrixes.y, sizeof(float) * (nnodes + 1)))
		throw DeviceMemoryAllocationException("could not allocate mtxd.y\n");	//CudaTest("couldn't allocate mtxd.y");
	if (cudaSuccess != cudaMalloc((void **)&devMatrixes.z, sizeof(float) * (nnodes + 1)))
		throw DeviceMemoryAllocationException("could not allocate mtxd.z\n");	//CudaTest("couldn't allocate mtxd.z");

	if (cudaSuccess != cudaMalloc((void **)&devMatrixes.x1, sizeof(float) * (nnodes + 1)))
		throw DeviceMemoryAllocationException("could not allocate mtxd.x1\n");	//CudaTest("couldn't allocate mtxd.x1");
	if (cudaSuccess != cudaMalloc((void **)&devMatrixes.y1, sizeof(float) * (nnodes + 1)))
		throw DeviceMemoryAllocationException("could not allocate mtxd.y1\n");	//CudaTest("couldn't allocate mtxd.y1");
	if (cudaSuccess != cudaMalloc((void **)&devMatrixes.z1, sizeof(float) * (nnodes + 1)))
		throw DeviceMemoryAllocationException("could not allocate mtxd.z1\n");	//CudaTest("couldn't allocate mtxd.z1");

	if (cudaSuccess != cudaMalloc((void **)&devMatrixes.theta, sizeof(float) * (nnodes + 1)))
		throw DeviceMemoryAllocationException("could not allocate mtxd.theta\n");	//CudaTest("couldn't allocate mtxd.theta");
	if (cudaSuccess != cudaMalloc((void **)&devMatrixes.phy, sizeof(float) * (nnodes + 1)))
		throw DeviceMemoryAllocationException("could not allocate mtxd.phy\n");	//CudaTest("couldn't allocate mtxd.phy");

	if (cudaSuccess != cudaMalloc((void **)&devMatrixes.theta1, sizeof(float) * (nnodes + 1)))
		throw DeviceMemoryAllocationException("could not allocate mtxd.theta1\n");	//CudaTest("couldn't allocate mtxd.theta1");
	if (cudaSuccess != cudaMalloc((void **)&devMatrixes.phy1, sizeof(float) * (nnodes + 1)))
		throw DeviceMemoryAllocationException("could not allocate mtxd.phy1\n");	//CudaTest("couldn't allocate mtxd.phy1");

	if (cudaSuccess != cudaMalloc((void **)&devMatrixes.ux, sizeof(float) * (nnodes + 1)))
		throw DeviceMemoryAllocationException("could not allocate mtxd.ux\n");	//CudaTest("couldn't allocate mtxd.ux");
	if (cudaSuccess != cudaMalloc((void **)&devMatrixes.uy, sizeof(float) * (nnodes + 1)))
		throw DeviceMemoryAllocationException("could not allocate mtxd.uy\n");	//CudaTest("couldn't allocate mtxd.uy");
	if (cudaSuccess != cudaMalloc((void **)&devMatrixes.uz, sizeof(float) * (nnodes + 1)))
		throw DeviceMemoryAllocationException("could not allocate mtxd.uz\n");	//CudaTest("couldn't allocate mtxd.uz");

	if (cudaSuccess != cudaMalloc((void **)&devMatrixes.deltaX, sizeof(float) * (nnodes + 1)))
		throw DeviceMemoryAllocationException("could not allocate mtxd.dx\n");	//CudaTest("couldn't allocate mtxd.dx");
	if (cudaSuccess != cudaMalloc((void **)&devMatrixes.deltaY, sizeof(float) * (nnodes + 1)))
		throw DeviceMemoryAllocationException("could not allocate mtxd.dy\n");	//CudaTest("couldn't allocate mtxd.dy");
	if (cudaSuccess != cudaMalloc((void **)&devMatrixes.deltaZ, sizeof(float) * (nnodes + 1)))
		throw DeviceMemoryAllocationException("could not allocate mtxd.dz\n");	//CudaTest("couldn't allocate mtxd.dz");

	if (cudaSuccess != cudaMalloc((void **)&devMatrixes.deltaTheta, sizeof(float) * (nnodes + 1)))
		throw DeviceMemoryAllocationException("could not allocate mtxd.dtheta\n");	//CudaTest("couldn't allocate mtxd.dtheta");
	if (cudaSuccess != cudaMalloc((void **)&devMatrixes.deltaPhy, sizeof(float) * (nnodes + 1)))
		throw DeviceMemoryAllocationException("could not allocate mtxd.dphy\n");	//CudaTest("couldn't allocate mtxd.dphy");

	if (cudaSuccess != cudaMalloc((void**)&devMatrixes.rndStates, sizeof(curandState) * blocks * FACTOR5 ))
			throw DeviceMemoryAllocationException("could not allocate rndStates\n");

	if (cudaSuccess != cudaMalloc((void **)&boxl.maxx, sizeof(float) * blocks * FACTOR1))
		throw DeviceMemoryAllocationException("could not allocate boxd.maxx\n");	//CudaTest("couldn't allocate boxd.maxx");
	if (cudaSuccess != cudaMalloc((void **)&boxl.maxy, sizeof(float) * blocks * FACTOR1))
		throw DeviceMemoryAllocationException("could not allocate boxd.maxy\n");	//CudaTest("couldn't allocate boxd.maxy");
	if (cudaSuccess != cudaMalloc((void **)&boxl.maxz, sizeof(float) * blocks * FACTOR1))
		throw DeviceMemoryAllocationException("could not allocate boxd.maxz\n");	//CudaTest("couldn't allocate boxd.maxz");
	if (cudaSuccess != cudaMalloc((void **)&boxl.minx, sizeof(float) * blocks * FACTOR1))
		throw DeviceMemoryAllocationException("could not allocate boxd.minx\n");	//CudaTest("couldn't allocate boxd.minx");
	if (cudaSuccess != cudaMalloc((void **)&boxl.miny, sizeof(float) * blocks * FACTOR1))
		throw DeviceMemoryAllocationException("could not allocate boxd.miny\n");	//CudaTest("couldn't allocate boxd.miny");
	if (cudaSuccess != cudaMalloc((void **)&boxl.minz, sizeof(float) * blocks * FACTOR1))
		throw DeviceMemoryAllocationException("could not allocate boxd.minz\n");	//CudaTest("couldn't allocate boxd.minz");

	std::cout << "Alloc matrixes CUDA finished " << std::endl;

	partPos.x = new float [cPar->nPart];
	partPos.y = new float [cPar->nPart];
	partPos.z = new float [cPar->nPart];
	partPos.theta = new float [ cPar->nPart];
	partPos.phy = new float [cPar->nPart];
	partPos.n = cPar->nPart;
}


void Fflcubh::initDevMatrixes() {
	std::cout << "Init matrixes" << std::endl;
	this->Fflcun2::initDevMatrixes();

	float *mass = new float [cPar->nPart];
	if (mass == NULL) throw HostMemoryAllocationException("cannot allocate mass");

	for(int i = 0; i < cPar->nPart; i++) {
		mass[i] = 1.0f;
	}
	if (cudaSuccess != cudaMemcpy(massl, mass, sizeof(float) * cPar->nPart, cudaMemcpyHostToDevice))
		throw DeviceMemoryCopyException("copying of mtx.mass to device failed\n");	//CudaTest("mtx.phy copy to device failed");
	delete[] mass;
}

void Fflcubh::freeMemory() {
	this->Fflcun2::freeMemory();

	cudaFree(massl);

	cudaFree(arrl.count);
	cudaFree(arrl.sort);
	cudaFree(arrl.start);
	cudaFree(arrl.err);
	cudaFree(arrl.child);

	cudaFree(boxl.maxx);
	cudaFree(boxl.maxy);
	cudaFree(boxl.maxz);
	cudaFree(boxl.minx);
	cudaFree(boxl.miny);
	cudaFree(boxl.minz);
}

Fflcun2::ParticlesPosition Fflcubh::run(int nSteps) {
	for(int j = 0; j < nSteps; j++) {

		if (cPar->ft == ConstParams::ROTATING) {
			chPar->hExtX = cPar->hExtXInit * sin(2 * M_PI * cPar->rff * chPar->time);
			chPar->hExtY = cPar->hExtYInit * cos(2 * M_PI * cPar->rff * chPar->time);
		}

		fillGloabalChangableBH(chPar);
		cudaThreadSynchronize();
#ifdef DEBUG
		cout << "Here "<<j<<endl;
#endif
		calcU(blocks);
		cudaThreadSynchronize();

#ifdef DEBUG
		cout << "Box "<<j<<endl;
#endif
		buildBox(blocks);
		cudaThreadSynchronize();

#ifdef DEBUG
		cout << "Tree "<<j<<endl;
#endif
		buildTree(blocks);
		cudaThreadSynchronize();

#ifdef DEBUG
		cout << "Summ "<<j<<endl;
#endif
		summarize(blocks);
		cudaThreadSynchronize();

#ifdef DEBUG
		cout << "Sort "<<j<<endl;
#endif
		sort(blocks);
		cudaThreadSynchronize();

#ifdef DEBUG
		cout << "Force "<<j<<endl;
#endif
		force(blocks);
		cudaThreadSynchronize();

#ifdef DEBUG
		cout << "Integrate "<<j<<endl;
#endif
		integrate(blocks);
		cudaThreadSynchronize();
		chPar->time += chPar->dTimeCurrent;
	}

	cudaMemcpy(partPos.x, devMatrixes.x, sizeof(float) * cPar->nPart, cudaMemcpyDeviceToHost);
	cudaMemcpy(partPos.y, devMatrixes.y, sizeof(float) * cPar->nPart, cudaMemcpyDeviceToHost);
	cudaMemcpy(partPos.z, devMatrixes.z, sizeof(float) * cPar->nPart, cudaMemcpyDeviceToHost);
	cudaMemcpy(partPos.theta, devMatrixes.theta, sizeof(float) * cPar->nPart, cudaMemcpyDeviceToHost);
	cudaMemcpy(partPos.phy, devMatrixes.phy, sizeof(float) * cPar->nPart, cudaMemcpyDeviceToHost);

	return partPos;
}

void Fflcubh::setN(int n) {
	this->freeMemory();
	cudaThreadSynchronize();

	cPar = new ConstParams;
	chPar = new ChangableParams;
	this->setSettings();
	this->cPar->nPart = n;

	nnodes = cPar->nPart * NODES_MULTY;
	if (nnodes < 1024 * blocks) nnodes = 1024 * blocks;
	while ((nnodes & (WARPSIZE - 1)) != 0) nnodes++;
	nnodes--;
	fprintf(stderr, "nnodes = %d, nbodies = %d\n", nnodes, cPar->nPart);

	this->allocMemory();
	this->initDevMatrixes();

	fillGloabalChangableBH(chPar);
	fillGloabalConstantBH(cPar);
	fillConstantPointers(cPar->nPart, nnodes, massl, arrl, boxl, devMatrixes);
	init(blocks, devMatrixes.rndStates);
	cudaThreadSynchronize();
}

float3 Fflcubh::getAverageU() {
	return averageU(blocks);
}

void Fflcubh::initRotStat() {
	if (devOldPhy != 0) this->destroyRotStat();
	if (cudaSuccess != cudaMalloc((void**)&devOldPhy, sizeof(float) * cPar->nPart ))
		throw DeviceMemoryAllocationException("Memory allocation of devOldPhy on device failed\n");

	getCurrentPhysBH(blocks, devOldPhy);
}


void Fflcubh::getRotStat(float* dPhy, int* pos) {
	getRingStatBH(blocks, nRings, cPar->nPart, &devOldPhy, dPhy, pos);
}


