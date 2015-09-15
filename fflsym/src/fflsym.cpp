//============================================================================
// Name        : fflsym.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <fflexprm.h>

using namespace std;

int main(int argc, char** argv) {
	time_t startTime, finishTime;
	std::string inputFileName;

	int rez=0;
	while ( (rez = getopt(argc,argv,"f:")) != -1){
		switch (rez){
		case 'f': inputFileName = optarg; break;
        }
	}
	if(inputFileName == "") inputFileName = "input.cfg";

	startTime = time(NULL);
	Experement e(inputFileName);
	e.initVisualization(/*Experement::CURRENT_POSITIONS |*/ Experement::STATISTICAL);
	e.run();

	finishTime = time(NULL);
	cout<<"Evaluation time is: "<< finishTime - startTime << "sec.";
	return 0;
}
