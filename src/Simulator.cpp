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
const string PARTICLES_PATH =  "C:/Users/Dallin/fluidsimstarter/data/particles/p.%03d";
const string GRIDS_PATH = "C:/Users/Dallin/fluidsimstarter/data/grids/grid.%03d";

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
			this->_grid_->applyParticleVelocities(this->_particles_);
			this->advectParticles(ts);
			this->_grid_->updateBuffer(this->_particles_, 1);
			this->_grid_->advectVelocity(ts);
			this->_grid_->applyExternalForces(ts, g);			
			this->_grid_->solvePressure(ts, FLUID_DENSITY, ATM_PRESSURE);
			this->_grid_->applyPressure(ts, FLUID_DENSITY);
			this->_grid_->extrapolateVelocity(1);
			this->_grid_->setSolidVelocities();
			// apply velocities to particles			
			this->applyGridVelsToP();
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
Eigen::Vector2d Simulator::interp(double x1, double xp, double x2, Eigen::Vector2d u1, Eigen::Vector2d u2) {
	int cell_size = double(this->_grid_->cellSize());
	// std::cout << "cell size: " << cell_size;
	double f1 = (xp-x1)/cell_size;
	double f2 = (x2-xp)/cell_size;
	// std::cout << "f values: " << f1 << ' ' << f2 << '\n';
	// std::cout << "u values: " << u1 << ' ' << u2 << '\n';
	// std::cout << "muliplied values: " << f2*u1 << ' ' << f1*u2 << '\n';

	return f1*u2 + f2*u1;
}

void Simulator::applyGridVelsToP() {
	for (Particle* p : this->_particles_) {
		Eigen::Vector2d newVel;
		this->_grid_->getVelocity(p->pos()[0], p->pos()[1], newVel);
		p->updateVel(newVel[0], newVel[1]); // incorrect for flip, but what the heck
	}
}

void Simulator::advectParticles(double t)
{
	Particle *p;
	Eigen::Vector2d velDiff, newPos, u;
	for(int i = 0; i < this->_particles_.size(); ++i)
	{
		p = this->_particles_[i];
		// this->_grid_->traceParticleDiff(p->pos()[0], p->pos()[1], t, velDiff);
		// p->updateVel(velDiff[0], velDiff[1]);

		this->_grid_->traceParticle(p->pos()[0], p->pos()[1], t, newPos);
		p->updatePos(newPos[0], newPos[1]);
		// find all grid cells in neighborhood
		// MacGrid* mcgriddle = this->_grid_;
		// int cell_size = mcgriddle->cellSize();
		// double px = p->pos()[0]/cell_size;
		// double py = p->pos()[1]/cell_size;
		// // std::cout << px << ' ' << py << '\n';
		// int gridx = int(px); // truncated
		// int gridy = int(py); // truncated
		// // std::cout << "u: " << gridx << " " << gridy << " " << mcgriddle->cellAt(gridx, gridy)->u() << '\n';

		// Eigen::Vector2d diff_v1 = mcgriddle->cellAt(gridx, gridy)->u() - mcgriddle->cellAt(gridx, gridy)->oldU();
		// Eigen::Vector2d diff_v2 = mcgriddle->cellAt(gridx+1, gridy)->u() - mcgriddle->cellAt(gridx+1, gridy)->oldU();
		// Eigen::Vector2d diff_v3 = mcgriddle->cellAt(gridx, gridy+1)->u() - mcgriddle->cellAt(gridx, gridy+1)->oldU();
		// Eigen::Vector2d diff_v4 = mcgriddle->cellAt(gridx+1, gridy+1)->u() - mcgriddle->cellAt(gridx+1, gridy+1)->oldU();

		// Eigen::Vector2d diff_interp1 = this->interp(gridx, px, gridx+1, diff_v1, diff_v2);
		// Eigen::Vector2d diff_interp2 = this->interp(gridx, px, gridx+1, diff_v3, diff_v4);
		// Eigen::Vector2d interp_diff_vel = this->interp(gridy, py, gridy+1, diff_interp1, diff_interp2);
		// // std::cout << "difference velocity: " << diff_v1 << '\n';
		// p->updateVel(interp_diff_vel[0]+0.5, interp_diff_vel[1]);
		// // std::cout << "pos: " << p->pos() << '\n';
		// std::cout << "vel: " << p->vel() << '\n';
		// p->updatePos(p->pos()[0]+p->vel()[0], p->pos()[1]+p->vel()[1]);
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
