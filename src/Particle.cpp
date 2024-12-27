/*
 * Particle.cpp
 *
 *  Created on: Jul 18, 2016
 *      Author: sflynn
 */

#include "Particle.h"

Particle::Particle(const double x, const double y) : _pos_(x, y)
{
    this->updateVel(0,0); // initialize velocity
    this->updateMass(20000); // nonzero mass by default
    this->updateVol(1);

    Eigen::Matrix2d I(2,2);
	I.setIdentity();
    this->updateDefGradP(I); // identity matrix in first frame
    this->updateDefGradE(I); // identity matrix in first frame


}

Particle::~Particle()
{

}

