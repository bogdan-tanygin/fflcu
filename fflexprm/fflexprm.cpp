/*
 * Experement.cpp
 *
 *  Created on: Jan 16, 2012
 *      Author: alexander
 */

#include "fflexprm.h"

Exprm::Task readTask(ConfigFile& conf)
{
	using namespace Exprm;
	string readTask;
	readTask = conf.read<string>("Task", "FOO_TASK");
	if (readTask == "KAPPA_TIME") return KAPPA_TIME;
	if (readTask == "SIMPLE_RUN") return SIMPLE_RUN;
	if (readTask == "RING_STAT") return RING_STAT;
	if (readTask == "HYST") return HYST;
	return FOO_TASK;	//todo - throw exception
}

Exprm::Variable readVariable(ConfigFile& conf)
{
	using namespace Exprm;
	string readVariable;
	readVariable = conf.read<string>("Variable", "FOO_VARIABLE");

	if (readVariable == "LX") return LX;
	if (readVariable == "LY") return LY;
	if (readVariable == "LZ") return LZ;

	if (readVariable == "HX") return HX;
	if (readVariable == "HY") return HY;
	if (readVariable == "HZ") return HZ;

	if (readVariable == "N_OF_PART") return N_OF_PART;
	if (readVariable == "ETA") return ETA;
	if (readVariable == "NYU") return NYU;

	return FOO_VARIABLE;
}

Exprm::Model readModel(ConfigFile& conf)
{
	using namespace Exprm;
	string readModel;
	readModel = conf.read<string>("Model", "FOO_MODEL");
	if (readModel == "N2") return N2;
	if (readModel == "BH") return BH;
	cout << "Warning: wrong model. N2 model has been set." << endl;
	return N2;
}


Experement::Experement(const std::string& confFileName) {
	inputFileName = confFileName;
	this->readConf();
	this->output.open(outputFileName.c_str(), std::fstream::app | std::fstream::out);
	this->visualization = 0;
}

Experement::~Experement() {
	//cout << "Experiment destructor" << endl;
	this->output.close();
	delete this->visualization;
}

void Experement::run() {
	const int nn = 100;
	long long globalNRun = 0;

	try{
		if (this->itsTask == Exprm::KAPPA_TIME) {

			for(; var <= varMax; var += varStep) {
				float3 kappaResult;
				kappaResult.x = 0;
				kappaResult.y = 0;
				kappaResult.z = 0;

				for (int it = 0; it < this->nIter; ++it) {
					float t = 0;
					float3 kappa;

					firstRun = true;

					if (mdl == Exprm::N2) sys = new Fflcun2(inputFileName);
					else sys = new Fflcubh(inputFileName);

					if (globalNRun == 0) {
						this->saveConf();
						sys->saveConf(this->output);
					}
					globalNRun++;

					this->setVar();
					firstRun = true;

					do {
						sys->run(nn);
						t += nn * dt;
						kappa = sys->getAverageU();

						if (this->visualization != 0) {
							if ((visGraphs & 0x01) == 0x01) {
								Fflcun2::ParticlesPosition pos;
								pos = sys->getState();
								visualization->showThis(pos.x, pos.y, pos.z, pos.theta, pos.phy, pos.n);
							}

							if ((visGraphs & 0x02) == 0x02) {
								//visualization->showCurrentState(t, it, var, kappa.x, kappa.y, kappa.z);
							}
						}

						if (t > tMax) {
							cout << "Not enough time" << endl;
							kappaResult.x = 0.0;
							kappaResult.y = 0.0;
							kappaResult.z = 0.0;
							break;
							//return;	//todo throw exception
						}
					} while (compareKappa(kappa) == false);

					kappaResult.x = (kappaResult.x * it + currentKappa.x) / (it + 1);
					kappaResult.y = (kappaResult.y * it + currentKappa.y) / (it + 1);
					kappaResult.z = (kappaResult.z * it + currentKappa.z) / (it + 1);
					delete sys;
					sys = 0;
				}

				output << var << "\t" << kappaResult.x << "\t"
						<< kappaResult.y << "\t" << kappaResult.z << std::endl;
			}
		}

		if (this->itsTask == Exprm::HYST) {
			int hystHalfLoopCntr = 0;
			float3 kappa;
			float t = 0.0;

			if (mdl == Exprm::N2) sys = new Fflcun2(inputFileName);
			else sys = new Fflcubh(inputFileName);

			this->saveConf();
			sys->saveConf(this->output);

			do {
				sys->run(nn);
				t += nn * dt;
				kappa = sys->getAverageU();

				if (this->visualization != 0) {
					if ((visGraphs & 0x01) == 0x01) {
						Fflcun2::ParticlesPosition pos;
						pos = sys->getState();
						visualization->showThis(pos.x, pos.y, pos.z, pos.theta, pos.phy, pos.n);
					}

					if ((visGraphs & 0x02) == 0x02) {
						//visualization->showCurrentState(t, it, var, kappa.x, kappa.y, kappa.z);
					}
				}

				if (t > tMax) {
					cout << "Not enough time" << endl;
					break;
				}
			} while (compareKappa(kappa) == false);


			while (hystHalfLoopCntr <= 2 * nHystLoops) {
				for (int it = 0; it < this->nIter; ++it) {
					float t = 0;
					t = 0.0f;

					this->setVar();
					cout << "var = " << var << endl;

					do {
						sys->run(nn);
						t += nn * dt;
						if (this->visualization != 0) {
							if ((visGraphs & 0x01) == 0x01) {
								Fflcun2::ParticlesPosition pos;
								pos = sys->getState();
								visualization->showThis(pos.x, pos.y, pos.z, pos.theta, pos.phy, pos.n);
							}

							if ((visGraphs & 0x02) == 0x02) {
								//visualization->showCurrentState(t, it, var, kappa.x, kappa.y, kappa.z);
							}
						}

					} while (t < tHyst);
					kappa = sys->getAverageU();

					output << var << "\t" << kappa.x << "\t"
							<< kappa.y << "\t" << kappa.z << std::endl;
				}

				var += varStep;
				if ((varStep > 0.0 && var > varMax) || (varStep < 0.0 && var < -varMax)) {
					varStep = varStep * (-1.0);
					hystHalfLoopCntr++;
				}
			}
		}

		if (this->itsTask == Exprm::SIMPLE_RUN) {
			float t = 0;
			time_t start, finish, evTime;
			float3 kappa;
			std::deque<double> averageVecMx, averageVecMy, averageVecMz;
			double averageMx = 0;
			double averageMy = 0;
			double averageMz = 0;

			if (mdl == Exprm::N2) sys = new Fflcun2(inputFileName);
			else sys = new Fflcubh(inputFileName);

			this->saveConf();
			sys->saveConf(this->output);

			while(t < tMax) {
				start = clock();
				sys->run(nn);
				finish = clock();
				evTime = ((finish - start) * 1000 )/ CLOCKS_PER_SEC;
				cout << "time for" << nn << "steps = " << evTime << "ms." << endl;

				t += nn * dt;

				if (this->visualization != 0) {
					if ((visGraphs & 0x01) == 0x01) {
						Fflcun2::ParticlesPosition pos;
						pos = sys->getState();
						visualization->showThis(pos.x, pos.y, pos.z, pos.theta, pos.phy, pos.n);
					}
					if ((visGraphs & 0x02) == 0x02) {
						float3 kappa = sys->getAverageU();//= getKappa(sys->getU());
						visualization->showCurrentState(t, 0, var, kappa.x, kappa.y, kappa.z);
					}
				}

				averageVecMx.push_back(kappa.x);
				averageVecMy.push_back(kappa.y);
				averageVecMz.push_back(kappa.z);

				if (averageVecMx.size() > averagingSteps) {
					averageMx = ((averageMx * (averageVecMx.size() - 1) - averageVecMx.front())
							+ averageVecMx.back()) / (averageVecMx.size() - 1);
					averageVecMx.pop_front();

					averageMy = ((averageMy * (averageVecMy.size() - 1) - averageVecMy.front())
							+ averageVecMy.back()) / (averageVecMy.size() -1);
					averageVecMy.pop_front();

					averageMz = ((averageMz * (averageVecMz.size() - 1) - averageVecMz.front())
							+ averageVecMz.back()) / (averageVecMz.size() - 1);
					averageVecMz.pop_front();
				} else {
					averageMx = (averageMx * (averageVecMx.size() - 1)
							+ averageVecMx.back()) / averageVecMx.size();
					averageMy = (averageMy * (averageVecMy.size() - 1)
							+ averageVecMy.back()) / averageVecMy.size();
					averageMz = (averageMz * (averageVecMz.size() - 1)
							+ averageVecMz.back()) / averageVecMz.size();
				}

				kappa = sys->getAverageU();
				output << t << "\t" << kappa.x << "\t" << kappa.y << "\t" << kappa.z
						<< "\t" << averageMx << "\t" << averageMy	<< "\t" << averageMz
						<< "\t" << averageVecMx.size() << std::endl;


			}
		}


		if (this->itsTask == Exprm::RING_STAT) {
			float t = 0;
			const int nrgs = 5;
			float* dphy = new float[nrgs];
			int* pd = new int[nrgs];
			int st = 0;
			const int dphyCheckStep = 100;

			if (mdl == Exprm::N2) sys = new Fflcun2(inputFileName);
			else sys = new Fflcubh(inputFileName);

			this->saveConf();
			sys->saveConf(this->output);
			sys->initRotStat();
			sys->setNRings(nrgs);


			while (t < tMax) {
				sys->run(nn);
				t += nn * dt;

				if (st % dphyCheckStep == 0) {
					sys->getRotStat(dphy, pd);
					output << t << "\t";
					for (int i = 0; i < nrgs; i++) {
						cout << i << "\t" << dphy[i] / (2 * M_PI * dt * nn * dphyCheckStep) << "\t" << pd[i] << endl;
					}
					for (int i = 0; i < nrgs; i++) {
						output << dphy[i] / (2 * M_PI * dt * nn * dphyCheckStep) << "\t";
					}
					for (int i = 0; i < nrgs; i++) {
						output << pd[i] << "\t";
					}
					output << endl;
					cout << "====================" << endl;
					st = 0;
				}
				st++;

				if (this->visualization != 0) {
					if ((visGraphs & 0x01) == 0x01) {
						Fflcun2::ParticlesPosition pos;
						pos = sys->getState();
						visualization->showThis(pos.x, pos.y, pos.z, pos.theta, pos.phy, pos.n);
					}
					if ((visGraphs & 0x02) == 0x02) {
						float3 kappa = sys->getAverageU();
						visualization->showCurrentState(t, 0, var, kappa.x, kappa.y, kappa.z);
					}
				}
			}
			delete[] dphy;
			delete[] pd;
			dphy = 0;
			pd = 0;
		}

		delete sys;
		sys = 0;
	} catch (DeviceMemoryException e) {
		cout << e.message << endl;
		return;
	}
}

void Experement::readConf() {
	ConfigFile conf(inputFileName);

	tMax = conf.read<float>("Tmax", 10.0);
	dt = conf.read<float>("dt", 0.01);
	nPart = conf.read<int>("N", 0);

	itsTask = readTask(conf);
	itsVar = readVariable(conf);

	outputFileName = conf.read<string>("FileName", "out.txt");
	nIter = conf.read<int>("Iterations", 1000);

	var = conf.read<float>("VariableMin", 3.0);
	varMax = conf.read<float>("VariableMax", 30.0);
	varStep = conf.read<float>("VariableStep", 1.0);
	kappaAcceptableInterval = conf.read<float>("KappaAcceptableInterval", 0.1);
	outputFileName = conf.read<std::string>("FileName", "output.txt");
	mdl = readModel(conf);
	tHyst = conf.read<float>("tHyst", 10.0);
	nHystLoops = conf.read<int>("nHystLoops", 1);
	averagingSteps = conf.read<int>("AveragingSteps", 10000);
}

bool Experement::compareKappa(float3 newKappa) {
	bool ret;
	if (firstRun == true) {
		currentKappa = newKappa;
		firstRun = false;
		return false;
	} else {
		if (fabs(newKappa.x - currentKappa.x) < kappaAcceptableInterval &&
			fabs(newKappa.y - currentKappa.y) < kappaAcceptableInterval &&
			fabs(newKappa.z - currentKappa.z) < kappaAcceptableInterval)
			ret = true;
		else ret = false;

		currentKappa = newKappa;
		return ret;
	}
}

void Experement::setVar() {
	switch (this->itsVar) {
		case Exprm::N_OF_PART:
			this->nPart = this->var;
			sys->setN(this->nPart);
			break;
		case Exprm::LX:
			sys->setLx(this->var);
			break;
		case Exprm::LY:
			sys->setLy(this->var);
			break;
		case Exprm::LZ:
			sys->setLz(this->var);
			break;
		case Exprm::HX:
			sys->setHx(this->var);
			break;
		case Exprm::HY:
			sys->setHy(this->var);
			break;
		case Exprm::HZ:
			sys->setHz(this->var);
			break;
		case Exprm::ETA:
			sys->setEta(this->var);
			break;
		case Exprm::NYU:
			sys->setNyu(this->var);
			break;
		default:
			break;
	}
}

void Experement::initVisualization(short vg) {
	this->visualization = new ConsoleVisualization();
	visGraphs = (short)vg;
}

void Experement::saveConf() {
	output << "# Task\t= " << (int)itsTask <<":";
	switch (itsTask) {
		case Exprm::KAPPA_TIME:
			output << "KAPPA_TIME" << endl;
			break;

		case Exprm::SIMPLE_RUN:
			output << "SIMPLE_RUN" << endl;
			break;

		default:
			output << "ERROR" << endl;
			break;
	}
	output << "# Var\t= " << (int)itsVar <<":";
	switch (itsVar) {
		case Exprm::ETA:
			output << "ETA" << endl;
			break;

		case Exprm::NYU:
			output << "NYU" << endl;
			break;

		case Exprm::HX:
			output << "HX" << endl;
			break;

		case Exprm::HY:
			output << "HY" << endl;
			break;

		case Exprm::HZ:
			output << "HZ" << endl;
			break;

		case Exprm::LX:
			output << "LX" << endl;
			break;

		case Exprm::LY:
			output << "LY" << endl;
			break;

		case Exprm::LZ:
			output << "LZ" << endl;
			break;

		case Exprm::N_OF_PART:
			output << "N_OF_PART" << endl;
			break;

		default:
			output << "ERROR" << endl;
			break;
	}
	output << "# VarMin\t= " << var << endl;
	output << "# VarMax\t= " << varMax << endl;
	output << "# nIter\t= " << nIter << endl;
	output << "# KappaInterval\t= " << kappaAcceptableInterval << endl;
	output << "# tMax\t= " << tMax << endl;
	output << "# AveragingSteps\t= " << averagingSteps << endl;
}
