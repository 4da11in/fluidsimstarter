//============================================================================
// Name        : fluidsimstarter.cpp
// Author      : Sean Flynn
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
using namespace std;

#include <eigen3/Eigen/Sparse>

#include "Simulator.h"
#include "MacGrid.h"
#include "GridCell.h"

const int NUM_FRAMES = 300;

int main()
{
	MacGrid *grid = new MacGrid(40, 30, 1);
	Simulator simulator(grid);
	simulator.run(NUM_FRAMES);
	return 0;
}
