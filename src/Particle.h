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
	void updatePos(const double x, const double y) {_pos_[0] = x; _pos_[1] = y;};

private:

	Eigen::Vector2d _pos_;
};

#endif /* PARTICLE_H_ */
