/*
 * Visualization.cpp
 *
 *  Created on: 12 сент. 2011
 *      Author: alexander
 */

#include "ConsoleVisualization.h"

ConsoleVisualization::ConsoleVisualization() {
	return;
}

/*virtual*/ void ConsoleVisualization::showThis(float* x, float* y, float* z, float* theta, float* phy, int n){
	cout<<"{";

	for(int i = 0; i < n; i++){
		cout<< setiosflags(ios::fixed) << setprecision(1) << x[i]<<";"<<y[i]<<";"<<z[i]<<";"<<theta[i]<<";"<<phy[i]<<";";
	}
	cout<<"}"<<endl;
}

/*virtual*/ void ConsoleVisualization::setSettings(int xScr, int yScr, float xSize, float ySize, float zSize){
	return;	//empty function
}

/*virtual*/ void ConsoleVisualization::showCurrentState(float t, int it, float var, float kx, float ky, float kz){
	cout << setiosflags(ios::fixed) << setprecision(2)
			<< "t:\t" << t <<
			"\nIteration:\t" << it <<
			"\nVariable:\t" << var <<
			"\nX-Susceptibility:\t" << kx <<
			"\nY-Susceptibility:\t" << ky <<
			"\nZ-Susceptibility:\t" << kz <<
			"\n=============================" << endl;
}

ConsoleVisualization::~ConsoleVisualization() {

}
