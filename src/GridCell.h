/*
 * GridCell.h
 *
 *  Created on: Jun 20, 2016
 *      Author: sflynn
 */

#ifndef GRIDCELL_H_
#define GRIDCELL_H_

#include <eigen3/Eigen/Sparse>
#include <iostream>
#include "Particle.h"

enum CellType {FLUID, AIR, SOLID, UNUSED};

class GridCell {
public:
	GridCell();
	virtual ~GridCell();

	//getters
	CellType const& type() const {return _type_;};
	int id() const {return _id_;};
	int const& layer() const {return _layer_;};
	Eigen::Vector2d& u() {return _u_;};
	Eigen::Vector2d& oldU() {return _oldU_;};
	Eigen::Vector2d& tempU() {return _tempU_;};
	double const& p() {return _p_;};

	std::vector<Particle*> getParticles() {return _particleNeighbors_;};


	//setters
	void setType(CellType type) {_type_ = type;};
	void setId(int id) {_id_ = id;};
	void setLayer(int layer) {_layer_ = layer;};
	void setU(Eigen::Vector2d u) {_u_[0] = u[0]; _u_[1] = u[1];};
	void setOldU(Eigen::Vector2d u) {_oldU_[0] = u[0]; _oldU_[1] = u[1];};
	void updateU(double x, double y) {_u_[0] = x; _u_[1] = y;};
	void setTempU(Eigen::Vector2d tempU) {_tempU_[0] = tempU[0]; _tempU_[1] = tempU[1];};
	void updateTempU(double x, double y) {_tempU_[0] = x; _tempU_[1] = y;};
	void setP(const double p) {_p_ = p;};
	void swapTempVelocity(void) {_u_[0] = _tempU_[0]; _u_[1] = _tempU_[1];};

	void addParticle(Particle* p){_particleNeighbors_.push_back(p);};


private:
	CellType _type_;

	int _id_;

	int _layer_;

	//velocity is stored at the center of the minimum cell edges
	Eigen::Vector2d _u_;
	Eigen::Vector2d _oldU_;

	Eigen::Vector2d _tempU_;

	//pressure is stored at the cell center
	double _p_;

	std::vector<Particle*> _particleNeighbors_;

};

#endif /* GRIDCELL_H_ */
