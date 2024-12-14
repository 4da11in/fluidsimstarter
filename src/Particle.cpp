/*
 * Particle.cpp
 *
 *  Created on: Jul 18, 2016
 *      Author: sflynn
 */

#include "Particle.h"

Particle::Particle(const double x, const double y) : _pos_(x, y)
{
    this->updateVel(0,0); // initialize velocity to zero
    this->updateMass(1); // nonzero mass by default
}

Particle::~Particle()
{

}

