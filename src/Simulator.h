/*
 * Simulator.h
 *
 *  Created on: Jul 6, 2016
 *      Author: sflynn
 */

#ifndef SIMULATOR_H_
#define SIMULATOR_H_

#include <iostream>
#include <vector>
#include <fstream>

#include "MacGrid.h"
#include "Particle.h"

#include <eigen3/Eigen/Sparse>

using namespace std;

class Simulator {
public:
	Simulator(MacGrid *grid);
	virtual ~Simulator();

	MacGrid* grid(void) const;
	vector<Particle*> particles(void) const;
	int const& fps() const {return _fps_;};

	void addParticle(Particle *p) {_particles_.push_back(p);};
	void setFps(const int fps) {_fps_ = fps;};

	void run(int frames);
	void addParticles(int count);
	void advectParticles(double t);
	void serializeGrids(int frame, const string path);
	void serializeParticles(int frame, const string path);

private:
	MacGrid *_grid_;
	vector<Particle*> _particles_;

	int _fps_;
};

#endif /* SIMULATOR_H_ */
