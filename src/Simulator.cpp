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
			// this->advectParticles(ts);
			this->pToGrid();
			if (curTime == 0) {
				this->initParticleVolumes();
			}
			this->_grid_->updateBuffer(this->_particles_, 1);
			this->_grid_->advectVelocity(ts);
			this->_grid_->applyExternalForces(ts, g);
			this->_grid_->solvePressure(ts, FLUID_DENSITY, ATM_PRESSURE);
			this->_grid_->applyPressure(ts, FLUID_DENSITY);
			this->_grid_->extrapolateVelocity(1);
			this->_grid_->setSolidVelocities();
			this->gridToP();
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
double h(double r) {
	if (r >= 0 && r <= 1) {
		return 1-r;
	} else if (r <= 0 && r >= -1) {
		return 1+r;
	} else {
		return 0;
	}
}
void Simulator::pToGrid()
{
	// in progress
	// MPM: apply mass and velocity to grid using proper interpolation
	
	for (int x = 0; x < this->_grid_->width(); ++x) {
		for (int y = 0; y < this->_grid_->height(); ++y) {
			GridCell* gc = this->_grid_->cellAt(x,y);
			Eigen::Vector2d vel(0,0);
			double cell_mass = 0;
			for (int i = 0; i < this->_particles_.size(); i++) {
				Particle* p = this->_particles_[i];
				// calculate weight using distance
				Eigen::Vector2d& particle_pos = p->pos();
				Eigen::Vector2d grid_pos(this->_grid_->cellSize()*x, this->_grid_->cellSize()*y);
				double weight = this->getWeight(x, y, particle_pos[0], particle_pos[1]);
				cell_mass += weight*p->mass();
				if (cell_mass != 0) {
					vel += p->vel()*p->mass()*weight/cell_mass;
				}
			}
			gc->setMass(cell_mass);
			gc->setOldU(vel);
			gc->setU(vel); // this will be updated by other steps later, oldU will remain the same
		}
	}
}

Eigen::Vector2d Simulator::interp(double x1, double xp, double x2, Eigen::Vector2d u1, Eigen::Vector2d u2) {
	int cell_size = double(this->_grid_->cellSize());
	double f1 = (xp-x1)/cell_size;
	double f2 = (x2-xp)/cell_size;

	return f1*u2 + f2*u1;
}

void Simulator::gridToP() {
	for (int x = 0; x < this->_grid_->width(); ++x) {
		for (int y = 0; y < this->_grid_->height(); ++y) { // loop through all cells
			GridCell* cell = this->_grid_->cellAt(x,y);
			Eigen::Vector2d diff = cell->u() - this->_grid_->cellAt(x,y)->oldU();
			cell->setDiff(diff);
		}
	}
	for (Particle* p : this->_particles_) {
		Eigen::Vector2d particle_pos = p->pos();
		int grid_cellx = particle_pos[0]/this->_grid_->cellSize();
		int grid_celly = particle_pos[1]/this->_grid_->cellSize();

		Eigen::Vector2d newVel = p->vel();
		
		for (int x = 0; x < this->_grid_->width(); ++x) {
			for (int y = 0; y < this->_grid_->height(); ++y) {
				double weight = this->getWeight(x, y, particle_pos[0], particle_pos[1]);
				Eigen::Vector2d diff = this->_grid_->cellAt(x,y)->diff();
				newVel += diff * weight;
			}
		}
		// zero velocity components if they are going into wall
		if (grid_cellx+1 > this->_grid_->width()) {
			newVel[0] = 0;	
		}
		if (grid_celly+1 > this->_grid_->height()) {
			newVel[1] = 0;
		}
		p->updateVel(newVel[0], newVel[1]);
	}
}

double Simulator::N(double x)
{
	x = abs(x);
	if (x >= 0 && x < 1) {
		return 0.5*pow(x,3) - pow(x,2) + 2.0/3.0;
	}
	if (x >= 1 && x < 2) {
		return -1.0/6.0 * pow(x,3) + pow(x,2) - 2*x + 4.0/3.0;
	}
    return 0.0;
}

double Simulator::getWeight(double i, double j, double x, double y)
{
    return N(x-i)*N(y-j);
}

void Simulator::advectParticles(double t)
{
	Particle *p;
	Eigen::Vector2d velDiff, newPos, u;
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

void Simulator::initParticleVolumes()
{
	for (int i = 0; i < this->_particles_.size(); i++) {
		Eigen::Vector2d p_pos = this->_particles_[i]->pos();
		double density = 0;
		for (int x = 0; x < this->_grid_->width(); ++x) {
			for (int y = 0; y < this->_grid_->height(); ++y) {
				double mass = this->_grid_->cellAt(x,y)->mass();
				double weight = this->getWeight(x, y, p_pos[0], p_pos[1]);
				density += mass * weight / pow(this->_grid_->cellSize(), 3);
			}
		}
		this->_particles_[i]->updateDens(density);
		this->_particles_[i]->updateVol(this->_particles_[i]->mass() / density);
	}
}
