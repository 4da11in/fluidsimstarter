/*
 * MacGrid.h
 *
 *  Created on: Jun 20, 2016
 *      Author: sflynn
 */

#ifndef UNIFORMGRID_H_
#define UNIFORMGRID_H_

#include <iostream>
#include <fstream>
#include <math.h>

#include "GridCell.h"
#include "Particle.h"

using namespace std;

class MacGrid {
public:
	MacGrid(int width, int height, int cellSize);
	virtual ~MacGrid();

	virtual double getMaxU(void);
	virtual double getMinCellSize(void);
	virtual int cellTypeCount(CellType type);
	virtual void updateBuffer(vector<Particle*> particles, int kcfl);
	virtual void setLayer(int layer);
	virtual void getVelocity(double x, double y, Eigen::Vector2d &result);
	virtual void advectVelocity(double t);
	virtual void traceParticle(double x, double y, double t, Eigen::Vector2d &result);
	virtual void swapTempVelocity(void);
	virtual void applyExternalForces(double t, double gravity);
	virtual void solvePressure(double t, double fluidDensity, double atmP);
	virtual void applyPressure(double t, double fluidDensity);
	virtual double getDivergence(int x, int y);
	virtual void extrapolateVelocity(int kcfl);
	virtual void setSolidVelocities(void);

	GridCell* cellAt(int i, int j) const;
	GridCell* cellAtWorldPos(double x, double y) const;
	void getNeighbors(int i, int j, GridCell** result) const;
	virtual int const& cellSize() const {return _cellSize_;};
	virtual int const& width() const {return _width_;};
	virtual int const& height() const {return _height_;};

	void serialize(const string path);

private:

	int _width_;
	int _height_;
	int _cellSize_;
	double _halfSize_;

	int _numCells_;

	GridCell ***_cells_;

	//pressure solve matrix coefficients
	Eigen::SparseMatrix<double>* _A_;

	//pressure vector, where the result of the pressure solve will go
	Eigen::VectorXd* _p_;

	//vector representing divergence of the velocity field
	Eigen::VectorXd* _b_;

	virtual double getInterpolatedValue(double x, double y, int index) const;
	virtual void getInterpWeights(double x, double y, int i, int j, double* result) const;
	virtual double getCellU(int i, int j, int index) const;
	virtual void getCellUComponents(int i, int j, int index, double* result) const;
	virtual int relabelFluidCells(void);
	virtual void buildPressureMatrix(double t, double fluidDensity, double atmP);

};

#endif /* UNIFORMGRID_H_ */
