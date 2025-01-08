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
const double PARTICLES_PER_FRAME = 300;//30;
const int EMIT_FRAMES = 1;//96;

//paths where serialized data will be written out (CHANGE THESE TO THE DESIRED PATH)
const string PARTICLES_PATH =  "C:/Users/Dallin/fluidsimstarter/data/particles/p.%03d";
const string GRIDS_PATH = "C:/Users/Dallin/fluidsimstarter/data/grids/grid.%03d";

Simulator::Simulator(MacGrid *grid) : _grid_(grid), _fps_(24)
{

}

Simulator::~Simulator()
{

}
double damping = 0;
double friction = 1;

void Simulator::run(int frames)
{

	double maxU, ts;
	double curTime = 0.0;
	double g = GRAVITY / this->_fps_;
	bool have_computed_volumes = false;

	bool printlots = false;
	printlots = true;

	for(int frame = 0; frame <= frames; ++frame)
	{
		cout << "Current frame: " << frame << "-----------------------------------------" << endl;
		if(frame < EMIT_FRAMES)
			this->addParticles(PARTICLES_PER_FRAME, frame);

		this->serializeGrids(frame, GRIDS_PATH);
		this->serializeParticles(frame, PARTICLES_PATH);

		while(curTime < frame)
		{
			if (printlots)
				cout << "\tcurrent time: " << curTime << endl;

			maxU = max(this->getMaxU(), this->_grid_->getMinCellSize());
			if (printlots)
				cout << "\t\tmaxU: " << maxU << endl;
			ts = this->_grid_->getMinCellSize() / maxU;
			ts = min(frame - curTime, ts);
			
			ts = min(ts, 1.0);
			ts = max(ts, 0.1);
			if (printlots)
				cout << "\t\ttimestep: " << ts << endl;
			
			this->pToGrid();
			Particle* p0 = this->_particles_[0];
			GridCell* mygc = this->_grid_->cellAt(7,22);
			if (!have_computed_volumes) {
				this->initParticleVolumes();
				have_computed_volumes = true;
			}
			// if (frame >= 30 && frame <= 60)
				// std::cout << "F used for external forces: \n" << p0->defGradE()*p0->defGradP() << '\n';
			this->_grid_->applyExternalForces(ts, g, this->_particles_, frame, damping, friction); // compute grid forces and update grid vels
			this->updateDeformationGradient(ts);
			this->gridToP(frame);
			this->advectParticles(ts);
			curTime += ts;

			if (printlots)
				cout << endl;
		}
	}
}

void Simulator::addParticles(int count, int frame)
{
	double randX, randY, randA, randR;
	Particle *p;
	for(int i = 0; i < count; ++i)
	{
		// 5-8, 21-24
		// randX = 5 + ((double)rand() / RAND_MAX) * 2;
		// randY = 21 + ((double)rand() / RAND_MAX) * 2;
		randA = ((double)rand() / RAND_MAX) * 2*3.14159;
		randR = sqrt(((double)rand() / RAND_MAX))*10;
		randX = 15 + randR*cos(randA);
		randY = 15 + randR*sin(randA);

		p = new Particle(randX, randY);
		this->addParticle(p);
		p->setParticleId(i + count*frame);
	}
}
double Simulator::getMaxU() {
	double maxu = 0;
	for (Particle* p : this->_particles_) {
		if (p->vel().norm() > maxu) {
			maxu = p->vel().norm();
		}
	}
	return maxu;
}
void Simulator::updateDeformationGradient(double t) {
	Eigen::Matrix2d I(2,2);
	I.setIdentity();
	for (int i = 0; i < this->_particles_.size(); i++) {
		Particle* p = this->_particles_[i];
		// if (i == 0)
			// std::cout << "\nparticle velocity before update def grad: \n" << p->vel() << '\n';
		Eigen::Matrix2d del_v = getDelV(p);
		Eigen::Matrix2d Fe = p->defGradE();
		Eigen::Matrix2d Fp = p->defGradP();
		Eigen::Matrix2d F_old = Fe*Fp;

		// std::cout << "det F: " << F_old.determinant() << ' ';
		
		// Eigen::Matrix2d dud_final{{1,0}, {0,10}};
		
		Eigen::Matrix2d F_final = (I + t*del_v)*F_old;//+ t*del_v
		// if (F_old(0,0) = 1 && F_old(1,1) == 1 && F_old(0,1) == 0 && F_old(1,0) == 0)
			// F_final = (I - t*del_v)*F_old;
			// F_final.transposeInPlace();
		// F_final(0,0) = 1; F_final(0,1) = 0; F_final(1,0) = 0; F_final(1,1) = 1.1;
		// F_final = dud_final;
		Eigen::Matrix2d Fe_temp = F_final*Fp.inverse();

		// compute svd
		Eigen::Matrix2d U, V, E;
		Eigen::JacobiSVD<Eigen::Matrix2d> svd(Fe_temp, Eigen::ComputeFullU | Eigen::ComputeFullV);

		U = svd.matrixU();
		Eigen::Vector2d singular_values = svd.singularValues();
		// clamp singular values
		// PARAMTERS
		double theta_c = 1.9e-1;
		double theta_s = 7.5e-2;
		for (double& d : singular_values) {
			d = std::clamp(d, 1-theta_c, 1+theta_s);
		}
		E = singular_values.asDiagonal();
		V = svd.matrixV();

		// Eigen::Matrix2d dud{{1,0}, {0,1}};
		// dud.setIdentity();
		Eigen::Matrix2d Fe_final = U*E*V.transpose();
		// Fe_final = dud;
		// Fe_final.setIdentity();
		p->updateDefGradE(Fe_final);
		
		// Eigen::Matrix2d dud2{{1.1,0}, {0,1}};
		Eigen::Matrix2d Fp_final = Fe_final.inverse()*F_final;
		// Fp_final.setIdentity();
		p->updateDefGradP(Fp_final);
		// if (p->getParticleId() == 0)
			// std::cout << "\n(t*del_v)\n" << (t*del_v) << '\n';
			// std::cout << "F: \n" << F_final << '\n';
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
			}
			gc->setMass(cell_mass);
			if (cell_mass != 0) {
				for (int i = 0; i < this->_particles_.size(); i++) {
					Particle* p = this->_particles_[i];
					Eigen::Vector2d& particle_pos = p->pos();
					double weight = this->getWeight(x, y, particle_pos[0], particle_pos[1]);
					vel += p->vel()*p->mass()*weight/cell_mass;
				}
			}
			// if (x == 7 && y == 22)
				// std::cout << "\nGrid vel after ptoGrid: \n" << vel;
			
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

void Simulator::gridToP(int frame) {
	for (int x = 0; x < this->_grid_->width(); ++x) {
		for (int y = 0; y < this->_grid_->height(); ++y) { // loop through all cells
			GridCell* cell = this->_grid_->cellAt(x,y);
			Eigen::Vector2d diff = cell->u() - cell->oldU();
			cell->setDiff(diff);
		}
	}
	for (Particle* p : this->_particles_) {
		// std::cout << p->deformation_gradient() << '\n';
		Eigen::Vector2d particle_pos = p->pos();
		int grid_cellx = particle_pos[0]/this->_grid_->cellSize();
		int grid_celly = particle_pos[1]/this->_grid_->cellSize();

		Eigen::Vector2d newVel = p->vel();
		
		for (int x = 0; x < this->_grid_->width(); ++x) {
			for (int y = 0; y < this->_grid_->height(); ++y) {
				double weight = this->getWeight(x, y, particle_pos[0], particle_pos[1]);				
				
				Eigen::Vector2d diff = this->_grid_->cellAt(x,y)->diff();
				newVel += diff * weight;
				if (diff[0] < 0) {
					// std::cout << '(' << x << ", " << y << ")  ";
				}

				// if (abs(weight) > 0.001 && p->getParticleId() == 0)
				// 	std::cout << "px, py: " << particle_pos[0] << ", " << particle_pos[1] << '\n';
					// std::cout << "x, y: " << x << ", " << y << '\n';
					// std::cout << "weight: " << weight << '\n';
					// std::cout << "diff: " << diff[0] << ", " << diff[1] << '\n';
				
			}
		}
		double maxvel = 3;
		if (newVel.norm() > maxvel)
		{
			newVel.normalize();
			newVel *= maxvel;
		}
		p->updateVel(newVel[0], newVel[1]);
	}
}

double Simulator::getWeight(double i, double j, double x, double y)
{
    return this->_grid_->N(x-i)*this->_grid_->N(y-j);
}

Eigen::Matrix2d Simulator::getDelV(Particle* p)
{
	double x = p->pos()[0];
	double y = p->pos()[1];
	Eigen::Matrix2d delv{{0,0},{0,0}};
	for (int i = 0; i < this->_grid_->width(); ++i) {
		for (int j = 0; j < this->_grid_->height(); ++j) {
			Eigen::Vector2d u = this->_grid_->cellAt(i, j)->u();
			Eigen::Vector2d dw = this->_grid_->getDelWeight(i, j, x, y);
			// if (p->getParticleId() == 0 && u.norm() > 0 && dw.norm() > 0)//i == 7 && j == 22 && 
				// std::cout << "cell: " << i << ' ' << j << ": " << "u: " << u.transpose() << " dw: " << dw.transpose() << '\n' << "\nu*dw^t \n" << u*dw.transpose() << '\n';
			// dw[1] = 0;
			// dw[0] = 0;
			delv += u*dw.transpose();
			// if (i == 7 && j == 22 && u.norm() > 0 && dw.norm() > 0)
				// std::cout << "particle " << p->getParticleId() << ": " << "\ndw:\n" << dw << "\n dw^t: \n" << dw.transpose() << '\n';
		}
	}

	// if (p->getParticleId() == 0)
		// std::cout << "\ndelv: \n" << delv;
	return delv;
}
void Simulator::advectParticles(double t)
{
	for(Particle *p : this->_particles_)
	{
		Eigen::Vector2d velDiff, newPos(0,0), u;
		// for (int i = 0; i < this->_grid_->width(); ++i) {
		// 	for (int j = 0; j < this->_grid_->height(); ++j) {
		// 		Eigen::Vector2d gridPos(i,j);
		// 		newPos += this->getWeight(i, j, p->pos()[0], p->pos()[1])*gridPos;
		// 	}
		// }
		newPos = p->pos() + t*p->vel();
		p->updatePos(newPos[0], newPos[1]);
	}
}

void Simulator::serializeGrids(int frame, const string path)
{
	char spath[300];
	sprintf(spath, path.c_str(), frame);
	// cout << "\twriting grid to " << spath << endl;
	this->_grid_->serialize(spath);
}

void Simulator::serializeParticles(int frame, const string path)
{
	char spath[300];
	sprintf(spath, path.c_str(), frame);
	// cout << "\twriting particles to " << spath << endl;
	ofstream f;
	f.open(spath);

	Eigen::Vector2d pos, u;
	for(int i = 0; i < this->_particles_.size(); ++i)
	{
		pos = this->_particles_[i]->pos();
		// this->_grid_->getVelocity(pos[0], pos[1], u);
		u = this->_particles_[i]->vel();
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
		// this->_particles_[i]->updateDens(density);
		this->_particles_[i]->updateVol(this->_particles_[i]->mass() / density);
	}
}
