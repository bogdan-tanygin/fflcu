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
 * Fflcun2.cpp
 *
 *  Created on: Jan 12, 2012
 *      Author: alexander
 */

#include "fflcu.h"
#define SQ(x) ((x)*(x))
const float myu0 = 4e-7 * M_PI;
const float kb = 1.3806488E-23;

inline float square(float x) {
	return x * x;
}

Fflcun2::Fflcun2(std::string fname): devOldPhy(0) {
	cPar = new ConstParams;
	chPar = new ChangableParams;
	configFileName = fname;
	this->setSettings();

	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, 0);
	//TODO - throw exception if no device found

	blocks = ceil(cPar->nPart * 1.0 / SHARED_ARRAY);
	blocks += blocks % deviceProp.multiProcessorCount;
	threads = cPar->nPart / blocks;
	cPar->nPart = threads * blocks;
	std::cout << "After correction number of particles = " << cPar->nPart << std::endl;

	srand(time(NULL));

	this->allocMemory();
	this->initDevMatrixes();

	fillGloabalChangable(chPar);
	fillGloabalConstant(cPar);

	runSetupKernel(blocks, threads, devMatrixes);
	cudaThreadSynchronize();
}

Fflcun2::Fflcun2(): devOldPhy(0) {
}

void Fflcun2::setConfig(std::string fname) {	//WARNING - test it
	this->freeMemory();
	cPar = new ConstParams;
	chPar = new ChangableParams;
	configFileName = fname;
	this->setSettings();

	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, 0);
	//TODO - throw exception if no device found

	blocks = ceil((float)cPar->nPart / SHARED_ARRAY);
	blocks += blocks % deviceProp.multiProcessorCount;
	threads = cPar->nPart / blocks;
	cPar->nPart = threads * blocks;
	std::cout << "After correction number of particles = " << cPar->nPart << std::endl;

	srand(time(NULL));

	this->allocMemory();
	this->initDevMatrixes();

	fillGloabalChangable(chPar);
	fillGloabalConstant(cPar);

	runSetupKernel(blocks, threads, devMatrixes);
	cudaThreadSynchronize();
}

Fflcun2::~Fflcun2() {
	std::cout << "N2 destructor" << std::endl;
	this->freeMemory();
}

void Fflcun2::allocMemory() {
	std::cout << "N2: Start allocation" << std::endl;

	cudaMalloc((void**)&devMatrixes.deltaPhy, sizeof(float) * cPar->nPart);
	cudaMalloc((void**)&devMatrixes.deltaTheta, sizeof(float) * cPar->nPart);
	cudaMalloc((void**)&devMatrixes.deltaX, sizeof(float) * cPar->nPart);
	cudaMalloc((void**)&devMatrixes.deltaY, sizeof(float) * cPar->nPart);
	cudaMalloc((void**)&devMatrixes.deltaZ, sizeof(float) * cPar->nPart);
	cudaMalloc((void**)&devMatrixes.phy, sizeof(float) * cPar->nPart);
	cudaMalloc((void**)&devMatrixes.theta, sizeof(float) * cPar->nPart);
	cudaMalloc((void**)&devMatrixes.phy1, sizeof(float) * cPar->nPart);
	cudaMalloc((void**)&devMatrixes.theta1, sizeof(float) * cPar->nPart);
	cudaMalloc((void**)&devMatrixes.ux, sizeof(float) * cPar->nPart);
	cudaMalloc((void**)&devMatrixes.uy, sizeof(float) * cPar->nPart);
	cudaMalloc((void**)&devMatrixes.uz, sizeof(float) * cPar->nPart);
	cudaMalloc((void**)&devMatrixes.x, sizeof(float) * cPar->nPart);
	cudaMalloc((void**)&devMatrixes.x1, sizeof(float) * cPar->nPart);
	cudaMalloc((void**)&devMatrixes.y, sizeof(float) * cPar->nPart);
	cudaMalloc((void**)&devMatrixes.y1, sizeof(float) * cPar->nPart);
	cudaMalloc((void**)&devMatrixes.z, sizeof(float) * cPar->nPart);
	cudaMalloc((void**)&devMatrixes.z1, sizeof(float) * cPar->nPart);
	cudaMalloc((void**)&devMatrixes.rndStates, sizeof(curandState) * cPar->nPart );

	std::cout << "N2: Device allocation succeeded" << std::endl;

	partPos.x = new float [cPar->nPart];
	if(partPos.x == NULL){
		throw HostMemoryAllocationException("partPos.x == NULL");
	}
	partPos.y = new float [cPar->nPart];
	if(partPos.y == NULL){
		throw HostMemoryAllocationException("partPos.y == NULL");
	}
	partPos.z = new float [cPar->nPart];
	if(partPos.z == NULL){
		throw HostMemoryAllocationException("partPos.z == NULL");
	}
	partPos.theta = new float [ cPar->nPart];
	if(partPos.theta == NULL){
		throw HostMemoryAllocationException("partPos.theta == NULL");
	}
	partPos.phy = new float [cPar->nPart];
	if(partPos.phy == NULL){
		throw HostMemoryAllocationException("partPos.phy == NULL");
	}
	partPos.n = cPar->nPart;

	std::cout << "N2: Host allocation succeeded" << std::endl;
}

void Fflcun2::freeMemory() {
	cudaFree(devMatrixes.deltaPhy);
	cudaFree(devMatrixes.deltaTheta);
	cudaFree(devMatrixes.deltaX);
	cudaFree(devMatrixes.deltaY);
	cudaFree(devMatrixes.deltaZ);
	cudaFree(devMatrixes.phy);
	cudaFree(devMatrixes.theta);
	cudaFree(devMatrixes.phy1);
	cudaFree(devMatrixes.theta1);
	cudaFree(devMatrixes.ux);
	cudaFree(devMatrixes.uy);
	cudaFree(devMatrixes.uz);
	cudaFree(devMatrixes.x);
	cudaFree(devMatrixes.x1);
	cudaFree(devMatrixes.y);
	cudaFree(devMatrixes.y1);
	cudaFree(devMatrixes.z);
	cudaFree(devMatrixes.z1);
	cudaFree(devMatrixes.rndStates);

	devMatrixes.deltaPhy = 0;
	devMatrixes.deltaTheta = 0;
	devMatrixes.deltaX = 0;
	devMatrixes.deltaY = 0;
	devMatrixes.deltaZ = 0;
	devMatrixes.phy = 0;
	devMatrixes.theta = 0;
	devMatrixes.phy1 = 0;
	devMatrixes.theta1 = 0;
	devMatrixes.ux = 0;
	devMatrixes.uy = 0;
	devMatrixes.uz = 0;
	devMatrixes.x = 0;
	devMatrixes.x1 = 0;
	devMatrixes.y = 0;
	devMatrixes.y1 = 0;
	devMatrixes.z = 0;
	devMatrixes.z1 = 0;
	devMatrixes.rndStates = 0;

	delete cPar;
	delete chPar;

	delete[] partPos.x;
	delete[] partPos.y;
	delete[] partPos.z;
	delete[] partPos.theta;
	delete[] partPos.phy;

	cPar = 0;
	chPar = 0;

	partPos.x = 0;
	partPos.y = 0;
	partPos.z = 0;
	partPos.theta = 0;
	partPos.phy = 0;

	this->destroyRotStat();
}

bool Fflcun2::checkParticlesLocation(unsigned int maxN) {
	if (cPar->cf == ConstParams::BARREL) {
		if (square(partPos.x[maxN] - cPar->barrelR) +
			square(partPos.y[maxN] - cPar->barrelR) > square(cPar->barrelR - 1.0f)) return false;
	}

	for(unsigned int i = 0; i < maxN; i++) {

		if (fabs(partPos.x[i] - partPos.x[maxN]) < 2.0f &&
			fabs(partPos.y[i] - partPos.y[maxN]) < 2.0f &&
			fabs(partPos.z[i] - partPos.z[maxN]) < 2.0f) {
			if (square(partPos.x[i] - partPos.x[maxN]) +
				square(partPos.y[i] - partPos.y[maxN]) +
				square(partPos.z[i] - partPos.z[maxN]) < 4.0f) {
				return false;
			}
		}
	}

	return true;
}

void Fflcun2::initDevMatrixes() {
	std::cout << "N2: Start filling" << std::endl;
	for (int i = 0; i < cPar->nPart; i++) {
		int currAttempt = 0;
		do {
			if (currAttempt >= this->attempts) {
				printf("Too dense system");
				exit(1);
			}
			currAttempt++;
			if (cPar->cf == ConstParams::BARREL) {
				partPos.x[i] = 1.0f + (cPar->barrelR * 2.0 - 2.0f) * rand() / RAND_MAX;
				partPos.y[i] = 1.0f + (cPar->barrelR * 2.0 - 2.0f) * rand() / RAND_MAX;
			} else {
				partPos.x[i] = 1.0f + (cPar->lx - 2.0f) * rand() / RAND_MAX;
				partPos.y[i] = 1.0f + (cPar->ly - 2.0f) * rand() / RAND_MAX;
			}

			partPos.z[i] = 1.0f + (cPar->lz - 2.0f) * rand() / RAND_MAX;
			//std::cout << "Attempt " << currAttempt << std::endl;
		} while(checkParticlesLocation(i) == false);
		//std::cout << "Set N " << i << " of " << cPar->nPart << std::endl;
	}

	for (int i = 0; i < cPar->nPart; i++) {
		partPos.theta[i] = M_PI * rand() / RAND_MAX;
		partPos.phy[i] = 2.0f * M_PI * rand() / RAND_MAX;
	}

	if (cudaSuccess != cudaMemcpy(devMatrixes.x, partPos.x, sizeof(float) * cPar->nPart, cudaMemcpyHostToDevice))
		throw DeviceMemoryCopyException("copying of mtx.x to device failed\n");	//CudaTest("mtx.x copy to device failed");
	if (cudaSuccess != cudaMemcpy(devMatrixes.y, partPos.y, sizeof(float) * cPar->nPart, cudaMemcpyHostToDevice))
		throw DeviceMemoryCopyException("copying of mtx.y to device failed\n");	//CudaTest("mtx.y copy to device failed");
	if (cudaSuccess != cudaMemcpy(devMatrixes.z, partPos.z, sizeof(float) * cPar->nPart, cudaMemcpyHostToDevice))
		throw DeviceMemoryCopyException("copying of mtx.z to device failed\n");	//CudaTest("mtx.z copy to device failed");
	if (cudaSuccess != cudaMemcpy(devMatrixes.theta, partPos.theta, sizeof(float) * cPar->nPart, cudaMemcpyHostToDevice))
		throw DeviceMemoryCopyException("copying of mtx.theta to device failed\n");	//CudaTest("mtx.theta copy to device failed");
	if (cudaSuccess != cudaMemcpy(devMatrixes.phy, partPos.phy, sizeof(float) * cPar->nPart, cudaMemcpyHostToDevice))
		throw DeviceMemoryCopyException("copying of mtx.phy to device failed\n");	//CudaTest("mtx.phy copy to device failed");

	if (cudaSuccess != cudaMemset(devMatrixes.x1, 0, sizeof(float) * cPar->nPart))
		throw DeviceMemsetException("setting of devMatrixes.x1 on device failed\n");	//CudaTest("devMatrixes.x1 set on device failed");
	if (cudaSuccess != cudaMemset(devMatrixes.y1, 0, sizeof(float) * cPar->nPart))
		throw DeviceMemsetException("setting of devMatrixes.y1 on device failed\n");	//CudaTest("devMatrixes.y1 set on device failed");
	if (cudaSuccess != cudaMemset(devMatrixes.z1, 0, sizeof(float) * cPar->nPart))
		throw DeviceMemsetException("setting of devMatrixes.z1 on device failed\n");	//CudaTest("devMatrixes.z1 set on device failed");
	if (cudaSuccess != cudaMemset(devMatrixes.theta1, 0, sizeof(float) * cPar->nPart))
		throw DeviceMemsetException("setting of devMatrixes.theta1 on device failed\n");	//CudaTest("devMatrixes.theta1 set on device failed");
	if (cudaSuccess != cudaMemset(devMatrixes.phy1, 0, sizeof(float) * cPar->nPart))
		throw DeviceMemsetException("setting of devMatrixes.phy1 on device failed\n");	//CudaTest("devMatrixes.phy1 set on device failed");
	std::cout << "N2: Filling succeeded" << std::endl;
}

void Fflcun2::setSettings() {
	ConfigFile conf(configFileName);

	float lx, ly, lz;	//x, y and z dimensions of the value
	float t2Ch, nyuAbs, etaAbs, v;
	int devices;

	lx = conf.read<double>("Lx", 2.2);
	ly = conf.read<double>("Ly", 60.0);
	lz = conf.read<double>("Lz", 60.0);

	if(lx <= 2.0 || ly <= 2.0 || lz <= 2.0) {
		std::cout << "Value size smaller then Particle";
		exit(1);
	}

	cPar->nPart = conf.read<int>("N", 0);

	cPar->lx = lx;
	cPar->ly = ly;
	cPar->lz = lz;
	cPar->timeMax = conf.read<double>("Tmax", 0.0);
/*	cPar->eta = conf.read<double>("eta", 0.1);	//trans
	cPar->nyu = conf.read<double>("nyu", 0.1);*/	//rot
	cPar->dTimeInitial = conf.read<double>("dt", 0.01);
	chPar->dTimeCurrent = cPar->dTimeInitial;
	chPar->sqrtdTime = sqrt(cPar->dTimeInitial);
	chPar->hExtX = conf.read<double>("hExtX", 0.0);
	chPar->hExtY = conf.read<double>("hExtY", 0.0);
	chPar->hExtZ = conf.read<double>("hExtZ", 0.0);
	chPar->time = 0.0;

	cPar->hExtXInit = chPar->hExtX;
	cPar->hExtYInit = chPar->hExtY;
	cPar->hExtZInit = chPar->hExtZ;

	cPar->gravitation = conf.read<bool>("gravitation", false);
	cPar->r = conf.read<double>("r", 1e-8);

	cPar->myu = conf.read<double>("myu", 1.4e5); //1.4*100000.0;
	cPar->roParticles = 8000;
	cPar->roEnvironment = 1000;
	cPar->thermalBath = conf.read<bool>("thermalBath", false);

	if(conf.read<string>("ContainerForm", "BRICK") == "BRICK")
		cPar->cf = ConstParams::BRICK;
	else cPar->cf = ConstParams::BARREL;

	cPar->barrelR = cPar->ly * 0.5f;
	if (cPar->cf == ConstParams::BARREL) cPar->lx = cPar->barrelR * 2.0f;

	this->attempts = conf.read<int>("Attempts", 1000);

	cPar->sigma = conf.read<double>("sigma", 2.1f);
	cPar->sigmaWall = conf.read<double>("sigmaWall", 1.0f);
	cPar->eps = conf.read<double>("eps", 0.0050f);
	cPar->epsWall = conf.read<double>("epsWall", 0.00010f);
	cPar->temper = conf.read<double>("T", 300.0f);

	cPar->rff = conf.read<double>("RotatingFieldFrequancy", 0.01);
	if(conf.read<string>("FieldType", "PERMANENT") == "PERMANENT")
		cPar->ft = ConstParams::PERMANENT;
	else cPar->ft = ConstParams::ROTATING;

	cPar->ljCutOffR = 3.0;
	cPar->ljCutOffForce = 24.0 * cPar->eps * (2.0 / pow(cPar->ljCutOffR, 13) - 1.0 / pow(cPar->ljCutOffR, 7)) / cPar->sigma;
	cPar->ljCutOffForceWall = 24.0 * cPar->epsWall * (2.0 / pow(cPar->ljCutOffR, 13) - 1.0 / pow(cPar->ljCutOffR, 7)) / cPar->sigmaWall;

	carrierViscosity = conf.read<double>("CarrierViscosity", 0.89e-3);

	etaAbs = 3.0 * M_PI * carrierViscosity * 2 * cPar->r;
	nyuAbs = M_PI * carrierViscosity * pow(2 * cPar->r, 3);
	std::cout << "Eta: "<< etaAbs << " Nyu: " << nyuAbs << std::endl;

	v = 4.0 / 3.0 * M_PI * powf(cPar->r, 3);
	t2Ch = (cPar->r * cPar->r * cPar->roParticles) /
				(4.0 / 3.0 * M_PI * myu0 * cPar->myu * cPar->myu);
	std::cout << "T2ch = " << t2Ch << std::endl;

	cPar->eta = etaAbs * sqrt(t2Ch) / (v * cPar->roParticles);	//trans
	cPar->nyu = nyuAbs * sqrt(t2Ch) / (v * cPar->r * cPar->r * cPar->roParticles);;	//rot
	std::cout << "Eta~: "<< cPar->eta << " Nyu~: " << cPar->nyu << std::endl;

	/*cPar->qt = t2Ch * sqrtf(2.0 * kb * cPar->temper * etaAbs) /
			(v * cPar->r * cPar->roParticles);
	cPar->qr = t2Ch * sqrtf(3.0 * kb * cPar->temper * nyuAbs) /
			(v * powf(cPar->r, 2) * cPar->roParticles);*/
	cPar->qt = t2Ch * sqrtf(2.0 * kb * cPar->temper * etaAbs / sqrtf(t2Ch)) /
			(v * cPar->r * cPar->roParticles);
	cPar->qr = t2Ch * sqrtf(3.0 * kb * cPar->temper * nyuAbs / sqrtf(t2Ch)) /
			(v * powf(cPar->r, 2) * cPar->roParticles);
	std::cout << "qt: "<< cPar->qt << " qr: " << cPar->qr << std::endl;


	/*t2Ch = (cPar->r * cPar->r * cPar->roParticles) /
			(4.0 / 3.0 * M_PI * myu0 * cPar->myu * cPar->myu);
	std::cout << "T2ch = " << t2Ch << std::endl;


	std::cout << nyuAbs << "  " << etaAbs <<"  ";
	cPar->qt = t2Ch * sqrtf(6.0 * kb * cPar->temper * etaAbs / sqrtf(t2Ch)) /
			(v * cPar->r * cPar->roParticles);
	cPar->qr = t2Ch * sqrtf(6.0 * kb * cPar->temper * nyuAbs / sqrtf(t2Ch)) /
			(v * powf(cPar->r, 2) * cPar->roParticles);
	std::cout << cPar->qt << "  " << cPar->qr <<"  ";
	*/

	/*t2Ch = (cPar->r * cPar->r * cPar->roParticles) /
			(4.0 / 3.0 * M_PI * myu0 * cPar->myu * cPar->myu);
	std::cout << "T2ch = " << t2Ch << std::endl;

	nyuAbs = cPar->nyu * v * cPar->r * cPar->r * cPar->roParticles / sqrtf(t2Ch);
	etaAbs = cPar->eta * v * cPar->roParticles / sqrtf(t2Ch);
	std::cout << nyuAbs << "  " << etaAbs <<"  ";
	cPar->qt = t2Ch * sqrtf(6.0 * kb * cPar->temper * etaAbs / sqrtf(t2Ch)) /
			(v * cPar->r * cPar->roParticles);
	cPar->qr = t2Ch * sqrtf(6.0 * kb * cPar->temper * nyuAbs / sqrtf(t2Ch)) /
			(v * powf(cPar->r, 2) * cPar->roParticles);
	std::cout << cPar->qt << "  " << cPar->qr <<"  ";
	*/

	deviceNubmer = conf.read<int>("CardNumber", 0);

	cudaGetDeviceCount(&devices);
	if (deviceNubmer + 1 > devices) {
		throw WrongDeviceNumberException(deviceNubmer);
	} else {
		cudaSetDevice(deviceNubmer);
	}
}

Fflcun2::ParticlesPosition Fflcun2::run(int nSteps) {
	for(int j = 0; j < nSteps; j++) {
		runCalcU(blocks, threads, devMatrixes);
		cudaThreadSynchronize();

		if (cPar->ft == ConstParams::ROTATING) {
			chPar->hExtX = cPar->hExtXInit * sin(2 * M_PI * cPar->rff * chPar->time);
			chPar->hExtY = cPar->hExtYInit * cos(2 * M_PI * cPar->rff * chPar->time);
		}
		fillGloabalChangable(chPar);
		cudaThreadSynchronize();

		runOneStep(blocks, threads, devMatrixes);
		cudaThreadSynchronize();

		runApplyDeltas(blocks, threads, devMatrixes);
		cudaThreadSynchronize();
		chPar->time += chPar->dTimeCurrent;
	}

/*	cudaMemcpy(partPos.x, devMatrixes.x, sizeof(float) * cPar->nPart, cudaMemcpyDeviceToHost);
	cudaMemcpy(partPos.y, devMatrixes.y, sizeof(float) * cPar->nPart, cudaMemcpyDeviceToHost);
	cudaMemcpy(partPos.z, devMatrixes.z, sizeof(float) * cPar->nPart, cudaMemcpyDeviceToHost);
	cudaMemcpy(partPos.theta, devMatrixes.theta, sizeof(float) * cPar->nPart, cudaMemcpyDeviceToHost);
	cudaMemcpy(partPos.phy, devMatrixes.phy, sizeof(float) * cPar->nPart, cudaMemcpyDeviceToHost);
*/
	cudaMemcpy(partPos.x, devMatrixes.x, sizeof(float) * cPar->nPart, cudaMemcpyDeviceToHost);
	cudaMemcpy(partPos.y, devMatrixes.y, sizeof(float) * cPar->nPart, cudaMemcpyDeviceToHost);
	cudaMemcpy(partPos.z, devMatrixes.z, sizeof(float) * cPar->nPart, cudaMemcpyDeviceToHost);
	cudaMemcpy(partPos.theta, devMatrixes.theta, sizeof(float) * cPar->nPart, cudaMemcpyDeviceToHost);
	cudaMemcpy(partPos.phy, devMatrixes.phy, sizeof(float) * cPar->nPart, cudaMemcpyDeviceToHost);

	return partPos;
}

void Fflcun2::setN(int n) {
	this->freeMemory();

	cPar = new ConstParams;
	chPar = new ChangableParams;
	this->setSettings();
	this->cPar->nPart = n;

	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, 0);

	blocks = ceil(cPar->nPart * 1.0 / SHARED_ARRAY);
	blocks += blocks % deviceProp.multiProcessorCount;
	threads = cPar->nPart / blocks;
	cPar->nPart = threads * blocks;
	std::cout << "After correction number of particles = " << cPar->nPart << std::endl;

	srand(time(NULL));

	this->allocMemory();
	this->initDevMatrixes();

	fillGloabalChangable(chPar);
	fillGloabalConstant(cPar);

	runSetupKernel(blocks, threads, devMatrixes);
	cudaThreadSynchronize();
}

void Fflcun2::saveConf(std::fstream& out) {
	out << "# N\t= " << cPar->nPart << std::endl;
	out << "# dt\t= " << cPar->dTimeInitial << std::endl;
	out << "# Lx\t= " << cPar->lx << std::endl;
	out << "# Ly\t= " << cPar->ly << std::endl;
	out << "# Lz\t= " << cPar->lz << std::endl;
	out << "# nyu\t= " << cPar->nyu << std::endl;
	out << "# eta\t= " << cPar->eta << std::endl;
	out << "# r\t= " << cPar->r << std::endl;
	out << "# myu\t= " << cPar->myu << std::endl;
	out << "# RoPart\t= " << cPar->roParticles<< std::endl;
	out << "# RoEnv\t= " << cPar->roEnvironment << std::endl;
	out << "# Gravitation\t= " << cPar->gravitation << std::endl;
	out << "# ThermalBath= " << cPar->thermalBath << std::endl;
	out << "# Hx\t= " << chPar->hExtX << std::endl;
	out << "# Hy\t= " << chPar->hExtY << std::endl;
	out << "# Hz\t= " << chPar->hExtZ << std::endl;
	out << "# T\t= " << cPar->temper << std::endl;

	if (cPar->cf == ConstParams::BARREL) out << "# ContainerForm\t= BARREL" << std::endl;
	else out << "# ContainerForm\t= BRICK" << std::endl;

	if (cPar->ft == ConstParams::ROTATING) out << "# FieldType\t= ROTATING" << std::endl;
	else out << "# FieldType\t= PERMANENT" << std::endl;
	out << "# RotatingFieldFrequency\t= " << cPar->rff << std::endl;

	out << "# Epsilon\t= " << cPar->eps << std::endl;
	out << "# EpsilonWall\t= " << cPar->epsWall << std::endl;
	out << "# Sigma\t= " << cPar->sigma << std::endl;
	out << "# SigmaWall\t= " << cPar->sigmaWall << std::endl;

	out << "# CarrierViscosity\t= " << carrierViscosity << std::endl;
}

float3 Fflcun2::getAverageU() {
	return averageU(blocks, threads, devMatrixes);
}

void Fflcun2::initRotStat() {
	if (devOldPhy != 0) this->destroyRotStat();
	if (cudaSuccess != cudaMalloc((void**)&devOldPhy, sizeof(float) * cPar->nPart ))
		throw DeviceMemoryAllocationException("Memory allocation of devOldPhy on device failed\n");

	getCurrentPhys(blocks, threads, devMatrixes, &devOldPhy);
}

void Fflcun2::destroyRotStat() {
	if (cudaSuccess !=cudaFree(devOldPhy))
		throw DeviceMemoryAllocationException("Free memory of devOldPhy on device failed\n");

	devOldPhy = 0;
}

void Fflcun2::getRotStat(float* dPhy, int* pos) {
	getRingStat(blocks, threads, devMatrixes, nRings, cPar->nPart, &devOldPhy, dPhy, pos);
}
