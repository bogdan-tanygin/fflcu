/*
 * Experement.h
 *
 *  Created on: Jan 16, 2012
 *      Author: alexander
 */

#ifndef FFLEXPRM_H_
#define FFLEXPRM_H_

#include <fstream>
#include <fflcu.h>
//#include <Fflcubh.h>
#include <vector_types.h>
#include <vector>
#include <list>
#include <deque>
#include <iterator>
#include "IVis.h"
#include "ConsoleVisualization.h"

namespace Exprm {
	enum Task{SIMPLE_RUN, KAPPA_TIME, RING_STAT, HYST, FOO_TASK};
	enum Variable{LX, LY, LZ, HX, HY, HZ, N_OF_PART, ETA, NYU, FOO_VARIABLE};
	enum Model{N2, BH, FOO_MODEL};
};

class Experement {
public:
	enum VisGraphsEnum{CURRENT_POSITIONS = 1, STATISTICAL = 2};	//1, 2, 4, 8...
	Experement(const std::string&);
	virtual ~Experement();
	void run();
	void initVisualization(short vg);
protected:
	IVis *visualization;
	short visGraphs;
private:
	void readConf();
	void saveConf();
	bool compareKappa(float3 newKappa);
	void setVar();

	Exprm::Task itsTask;
	Exprm::Variable itsVar;
	Exprm::Model mdl;
	float var;
	float varMax;
	float varStep;
	int nIter;
	float kappaAcceptableInterval;
	float tMax;
	float dt;
	float tHyst;	//time for one hyst step
	int nPart;
	int nHystLoops;
	float3 currentKappa;
	unsigned long averagingSteps;

	bool firstRun;

	std::string inputFileName;
	std::string outputFileName;
	std::fstream output;
	Fflcun2::ParticlesPosition currState;
	Fflcun2 *sys;

};

#endif /* FFLEXPRM_H_ */
