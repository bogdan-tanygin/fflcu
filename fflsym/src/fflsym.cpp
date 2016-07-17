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
