/*
 * Particle.h
 *
 *  Created on: Jul 18, 2016
 *      Author: sflynn
 */

#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/SVD>

class Particle {
public:
	Particle(const double x, const double y);
	virtual ~Particle();

	Eigen::Vector2d& pos(void) {return _pos_;};
	Eigen::Vector2d& vel(void) {return _vel_;};
	double mass(void) {return _mass_;};
	double dens(void) {return _density_;};
	double vol(void) {return _volume_;};
	int getParticleId() {return _id_;};

	Eigen::Matrix2d& defGradE(void) {return _deformation_gradient_e_;};
	Eigen::Matrix2d& defGradP(void) {return _deformation_gradient_p_;};


	void updatePos(const double x, const double y) {_pos_[0] = x; _pos_[1] = y;};
	void updateVel(const double x, const double y) {_vel_[0] = x; _vel_[1] = y;};
	
	void updateMass(const double m) {_mass_ = m;};
	void updateDens(const double d) {_density_ = d;};
	void updateVol(const double v) {_volume_ = v;};

	void updateDefGradP(Eigen::Matrix2d m) {_deformation_gradient_p_ = m;};	
	void updateDefGradE(Eigen::Matrix2d m) {_deformation_gradient_e_ = m;};	

	void setParticleId(int id) {_id_ = id;};
private:

	Eigen::Vector2d _pos_;
	Eigen::Vector2d _vel_;

	Eigen::Matrix2d _deformation_gradient_p_;
	Eigen::Matrix2d _deformation_gradient_e_;

	double _mass_;
	double _density_;
	double _volume_;

	int _id_;

};

#endif /* PARTICLE_H_ */
