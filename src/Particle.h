/*
 * Particle.h
 *
 *  Created on: Jul 18, 2016
 *      Author: sflynn
 */

#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <eigen3/Eigen/Sparse>

class Particle {
public:
	Particle(const double x, const double y);
	virtual ~Particle();

	Eigen::Vector2d& pos(void) {return _pos_;};
	Eigen::Vector2d& vel(void) {return _vel_;};
	double mass(void) {return _mass_;};
	double dens(void) {return _mass_;};
	double vol(void) {return _mass_;};
	Eigen::Vector2d& deformation_gradient(void) {return _vel_;};

	void updatePos(const double x, const double y) {_pos_[0] = x; _pos_[1] = y;};
	void updateVel(const double x, const double y) {_vel_[0] = x; _vel_[1] = y;};
	
	void updateMass(const double m) {_mass_ = m;};
	void updateDens(const double d) {_density_ = d;};
	void updateVol(const double v) {_volume_ = v;};
	void updateDefGrad(Eigen::Vector2d g) {_deformation_gradient_ = g;};

private:

	Eigen::Vector2d _pos_;
	Eigen::Vector2d _vel_;
	Eigen::Vector2d _deformation_gradient_;

	double _mass_;
	double _density_;
	double _volume_;

};

#endif /* PARTICLE_H_ */
