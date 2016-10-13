/*
 * Simulator.cpp
 *
 *  Created on: Jul 6, 2016
 *      Author: sflynn
 */

#include "Simulator.h"

//gravity constant, one grid cell is one meter square
const double GRAVITY = -9.8321849378;

const double FLUID_DENSITY = 1.0;
const double ATM_PRESSURE = 1.0;

//constants controlling how particles are added to the simulation
const double PARTICLES_PER_FRAME = 30;
const int EMIT_FRAMES = 96;

//paths where serialized data will be written out (CHANGE THESE TO THE DESIRED PATH)
const string PARTICLES_PATH =  "<YOUR_PATH_HERE>/particles/p.%03d";
const string GRIDS_PATH = "<YOUR_PATH_HERE>/grids/grid.%03d";

Simulator::Simulator(MacGrid *grid) : _grid_(grid), _fps_(24)
{

}

Simulator::~Simulator()
{

}

void Simulator::run(int frames)
{

	double maxU, ts;
	double curTime = 0.0;
	double g = GRAVITY / this->_fps_;

	for(int frame = 0; frame <= frames; ++frame)
	{
		cout << "Current frame: " << frame << "-----------------------------------------" << endl;
		if(frame < EMIT_FRAMES)
			this->addParticles(PARTICLES_PER_FRAME);

		this->serializeGrids(frame, GRIDS_PATH);
		this->serializeParticles(frame, PARTICLES_PATH);

		while(curTime < frame)
		{
			cout << "\tcurrent time: " << curTime << endl;

			maxU = max(this->_grid_->getMaxU(), this->_grid_->getMinCellSize());
			cout << "\t\tmaxU: " << maxU << endl;
			ts = this->_grid_->getMinCellSize() / maxU;
			ts = min(frame - curTime, ts);
			cout << "\t\ttimestep: " << ts << endl;
			this->advectParticles(ts);

			this->_grid_->updateBuffer(this->_particles_, 1);
			this->_grid_->advectVelocity(ts);
			this->_grid_->applyExternalForces(ts, g);
			this->_grid_->solvePressure(ts, FLUID_DENSITY, ATM_PRESSURE);
			this->_grid_->applyPressure(ts, FLUID_DENSITY);
			this->_grid_->extrapolateVelocity(1);
			this->_grid_->setSolidVelocities();

			curTime += ts;

			cout << endl;
		}
	}
}

void Simulator::addParticles(int count)
{
	double randX, randY;
	Particle *p;
	for(int i = 0; i < count; ++i)
	{
		randX = ((this->_grid_->width() / 1.7) * this->_grid_->cellSize() -
					   ((double)rand() / RAND_MAX) * 2.5);
		randY = ((this->_grid_->height() / 1.2) * this->_grid_->cellSize() +
				       ((double)rand() / RAND_MAX) * 1.5);
		p = new Particle(randX, randY);
		this->addParticle(p);
	}
}

void Simulator::advectParticles(double t)
{
	Particle *p;
	Eigen::Vector2d newPos, u;
	for(int i = 0; i < this->_particles_.size(); ++i)
	{
		p = this->_particles_[i];
		this->_grid_->traceParticle(p->pos()[0], p->pos()[1], t, newPos);
		p->updatePos(newPos[0], newPos[1]);
	}
}

void Simulator::serializeGrids(int frame, const string path)
{
	char spath[300];
	sprintf(spath, path.c_str(), frame);
	cout << "\twriting grid to " << spath << endl;
	this->_grid_->serialize(spath);
}

void Simulator::serializeParticles(int frame, const string path)
{
	char spath[300];
	sprintf(spath, path.c_str(), frame);
	cout << "\twriting particles to " << spath << endl;
	ofstream f;
	f.open(spath);

	Eigen::Vector2d pos, u;
	for(int i = 0; i < this->_particles_.size(); ++i)
	{
		pos = this->_particles_[i]->pos();
		this->_grid_->getVelocity(pos[0], pos[1], u);
		f << pos[0] << " " << pos[1] << " " << u[0] << " " << u[1] << endl;
	}
}
