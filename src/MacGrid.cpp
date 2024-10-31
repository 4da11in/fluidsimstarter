/*
 * MacGrid.cpp
 *
 *  Created on: Jun 20, 2016
 *      Author: sflynn
 */

#include "MacGrid.h"

const double NON_EXISTENT_VEL = 0.0;

MacGrid::MacGrid(int width, int height, int cellSize) :
	_width_(width), _height_(height), _cellSize_(cellSize)
{
	this->_numCells_ = width * height;
	this->_cells_ = new GridCell**[width];
	this->_halfSize_ = cellSize / 2.0;

	//initialize _cells_ to empty GridCells
	for(int x = 0; x < width; ++x)
	{
		this->_cells_[x] = new GridCell*[height];
		for(int y = 0; y < height; ++y)
		{
			this->_cells_[x][y] = new GridCell();
		}
	}

	//initialize the data structures for the pressure solve
	this->_A_ = new Eigen::SparseMatrix<double>(this->_numCells_, this->_numCells_);
	this->_p_ = new Eigen::VectorXd(this->_numCells_);
	this->_b_ = new Eigen::VectorXd(this->_numCells_);
}

MacGrid::~MacGrid()
{

}


double MacGrid::getMaxU(void)
{
	//get the maximum velocity component in the grid
	double maxU = 0.0;
	double mag = 0.0;
	for(int x = 0; x < this->_width_; ++x)
	{
		for(int y = 0; y < this->_height_; ++y)
		{
			mag = this->_cells_[x][y]->u().norm();
			if(mag > maxU)
				maxU = mag;
		}
	}
	return maxU;
}

double MacGrid::getMinCellSize(void)
{
	return this->_cellSize_;
}

int MacGrid::cellTypeCount(CellType type)
{
	//returns the number of cells of type type
	int count = 0;
	for(int x = 0; x < this->_width_; ++x)
		for(int y = 0; y < this->_height_; ++y)
			if(this->cellAt(x, y)->type() == type)
				++count;

	return count;
}

void MacGrid::updateBuffer(vector<Particle*> particles, int kcfl)
{
	//sets the cells that have fluid particles in them to type fluid, and creates a buffer of air
	//cells around the fluid cells

	this->setLayer(-1);

	Eigen::Vector2d p;
	GridCell *cell, *neighbor;
	GridCell *neighbors[4];

	//set cells with fluid
	const int numParticles = particles.size();
	for(int i = 0; i < numParticles; ++i)
	{
		p = particles[i]->pos();
		cell = this->cellAtWorldPos(p[0], p[1]);
		if(cell != NULL && cell->type() != SOLID)
		{
			cell->setType(FLUID);
			cell->setLayer(0);
		}
	}

	int maxLayer = max(2, kcfl);

	CellType type;
	int layer;

	//create a buffer of air around the fluid cells
	for(int index = 1; index < maxLayer; ++index)
	{
		for(int x = 0; x < this->_width_; ++x)
		{
			for(int y = 0; y < this->_height_; ++y)
			{
				cell = this->cellAt(x, y);
				type = cell->type();
				layer = cell->layer();
				if(layer == index - 1 && (type == FLUID || type == AIR))
				{
					this->getNeighbors(x, y, neighbors);

					for(int n = 0; n < 4; ++n)
					{
						neighbor = neighbors[n];
						if(neighbor != NULL && neighbor->layer() == -1 &&
						   neighbor->type() != SOLID)
						{
							neighbor->setType(AIR);
							neighbor->setLayer(index);
						}

					}
				}
			}
		}
	}

	for(int x = 0; x < this->_width_; ++x)
	{
		for(int y = 0; y < this->_height_; ++y)
		{
			cell = this->cellAt(x, y);
			if(cell->type() != SOLID && cell->layer() == -1)
			{
				cell->setType(UNUSED);
			}
		}
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

void MacGrid::applyParticleVelocities(vector<Particle *> particles)
{
	// in progress
	// need to apply velocities from particles to grid.
	// for each particle, 

	// for (int pnum = 0; pnum < particles.size(); ++pnum) {
	// 	// add this particle to the "particle neighbor" list of every grid cell that may be affected by this particle's velocity
	// 	Particle* p = particles[pnum];
	// 	double x = p->pos()[0];
	// 	double y = p->pos()[1];
	// 	int gridx = x/this->cellSize();
	// 	int gridy = y/this->cellSize();
	// 	// which cells might this particle affect?
	// 	// anything within 1 grid cell on the x (keeping in mind we may be at either extreme, and are rounding left and down)
	// 	this->cellAt(gridx, gridy)->addParticle(p);
	// 	if (gridx+1 < width())
	// 		this->cellAt(gridx+1, gridy)->addParticle(p);
	// 	if (gridx > 0)
	// 		this->cellAt(gridx-1, gridy)->addParticle(p);
	// 	// checking on the y
	// 	if (gridy+1 < height())
	// 		this->cellAt(gridx, gridy+1)->addParticle(p);
	// 	if (gridy > 0)
	// 		this->cellAt(gridx, gridy-1)->addParticle(p);
	// }

	for (int x = 0; x < width(); ++x) {
		for (int y = 0; y < height(); ++y) {
			GridCell* gc = this->cellAt(x,y);
			double W = 0;
			Eigen::Vector2d vel(0,0);
			for (int i = 0; i < particles.size(); i++) { // gc->getParticles()
			// should be computing per component
				Particle* p = particles[i];
				// calculate weight using distance (distance reduces weight, something like 1-d)
				Eigen::Vector2d& particle_pos = p->pos();
				Eigen::Vector2d grid_pos(cellSize()*x, cellSize()*y);
				// check if close enough in the x
				// need to shift assuming that grid positions are in center, not lower left
				// double x_component_y_pos = grid_pos[1]+0.5*double(cellSize());
				// bool x_close_enough = abs(particle_pos[0]-grid_pos[0]) < cellSize() && abs(particle_pos[1]-x_component_y_pos) < cellSize();
				// double y_component_x_pos = grid_pos[0]+0.5*double(cellSize());
				// double y_close_enough = abs(particle_pos[1]-grid_pos[1]) < cellSize() && abs(particle_pos[0]-y_component_x_pos) < cellSize();
				double w = 0;
				// if (x_close_enough && y_close_enough) {
				w = h((particle_pos[0]-grid_pos[0])/cellSize())*h((particle_pos[1]-grid_pos[1])/cellSize());
				// }
				W += w;
				vel += w*p->vel();
			}
			// std::cout << "W (total) is: " << W << '\n';
			// std::cout << x << " " << y<< " " << " " << vel<< " " << '\n';
			// use weights on velocities, sum all and divide by sum of weights. save this new velocity onto the grid.
			if (W != 0) {
				vel/=W;
			}
			gc->setOldU(vel);
			gc->setU(vel); // this will be updated by other steps later, oldU will remain the same
		}
	}
}

void MacGrid::setLayer(int layer)
{
	//sets the layer of all the cells in the grid to layer
	for(int x = 0; x < this->_width_; ++x)
	{
		for(int y = 0; y < this->_height_; ++y)
		{
			this->_cells_[x][y]->setLayer(layer);
		}
	}
}

void MacGrid::getVelocity(double x, double y, Eigen::Vector2d &result)
{
	//x and y are in world space, not grid space

	double dx = x / this->_cellSize_;
	double dy = y / this->_cellSize_;

	double hdx = dx - 0.5;
	double hdy = dy - 0.5;

	result[0] = this->getInterpolatedValue(dx, hdy, 0);
	result[1] = this->getInterpolatedValue(hdx, dy, 1);
}

void MacGrid::getVelocityDiff(double x, double y, Eigen::Vector2d &result)
{
	//x and y are in world space, not grid space

	double dx = x / this->_cellSize_;
	double dy = y / this->_cellSize_;

	double hdx = dx - 0.5;
	double hdy = dy - 0.5;

	result[0] = this->getInterpolatedDiff(dx, hdy, 0);
	result[1] = this->getInterpolatedDiff(hdx, dy, 1);
}

void MacGrid::advectVelocity(double t)
{
	//advects the velocity field using the backward particle trace

	double nt = t * -1.0;
	GridCell *cell;

	Eigen::Vector2d prevX, prevY, prevUX, prevUY;

	for(int x = 0; x < this->_width_; ++x)
	{
		for(int y = 0; y < this->_height_; ++y)
		{
			cell = this->cellAt(x, y);

			this->traceParticle(x * this->_cellSize_, y * this->_cellSize_ + this->_halfSize_,
							    nt, prevX);
			this->traceParticle(x * this->_cellSize_ + this->_halfSize_, y * this->_cellSize_,
								nt, prevY);

			this->getVelocity(prevX[0], prevX[1], prevUX);
			this->getVelocity(prevY[0], prevY[1], prevUY);

			cell->updateTempU(prevUX[0], prevUY[1]);
			// cell->setOldU(cell->u());
		}
	}

	this->swapTempVelocity();
}

void MacGrid::traceParticle(double x, double y, double t, Eigen::Vector2d &result)
{

	//x and y are in world space, not grid space

	Eigen::Vector2d v;

	//store velocity at current location in v
	this->getVelocity(x, y, v);

	//double ht = this->_halfSize_ * t;
	double ht = 0.5 * t;

	//store velocity at location half a timestep ago in v
	this->getVelocity(x + ht * v[0], y + ht * v[1], v);


	//advect by the velocity half a timestep ago
	result[0] = x + t * v[0];
	result[1] = y + t * v[1];
}

void MacGrid::traceParticleDiff(double x, double y, double t, Eigen::Vector2d &result)
{

	//x and y are in world space, not grid space

	Eigen::Vector2d v;

	//store velocity at current location in v
	this->getVelocityDiff(x, y, v);

	//double ht = this->_halfSize_ * t;
	double ht = 0.5 * t;

	//store velocity at location half a timestep ago in v
	this->getVelocityDiff(x + ht * v[0], y + ht * v[1], v);


	//advect by the velocity half a timestep ago
	result[0] = x + t * v[0];
	result[1] = y + t * v[1];
}

void MacGrid::swapTempVelocity(void)
{
	for(int x = 0; x < this->_width_; ++x)
	{
		for(int y = 0; y < this->_height_; ++y)
		{
			this->cellAt(x, y)->swapTempVelocity();
		}
	}
}

void MacGrid::applyExternalForces(double t, double gravity)
{
	GridCell *cell, *n;
	double vg = gravity * t;

	for(int x = 0; x < this->_width_; ++x)
	{
		for(int y = 0; y < this->_height_; ++y)
		{
			cell = this->cellAt(x, y);

			//update the velocity if it borders a FLUID cell
			if(cell->type() == FLUID)
			{
				cell->u()[1] += vg;
			}
			else
			{
				n = this->cellAt(x, y - 1);
				if(n != NULL and n->type() == FLUID)
					cell->u()[1] += vg;
			}
		}
	}
}

void MacGrid::solvePressure(double t, double fluidDensity, double atmP)
{
	this->buildPressureMatrix(t, fluidDensity, atmP);

	//use conjugate gradient for the matrix solve
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > cg;
	cg.compute(*this->_A_);
	*this->_p_ = cg.solve(*this->_b_);
}

void MacGrid::applyPressure(double t, double fluidDensity)
{
	cout << "applyPressure: DONE!!!" << endl;
		
	for (int x = 0; x < width(); ++x) {
		for (int y = 0; y < height(); ++y) { // loop through all cells
			GridCell* cell = this->cellAt(x, y);
			Eigen::VectorXd p = *(this->_p_);
			// x gradient
			float x_grad = 0;
			if (cellAt(x-1, y) != NULL) {
				bool borders_fluid_cell = cellAt(x-1, y)->type()==FLUID || cellAt(x, y)->type()==FLUID;
				bool borders_solid_cell = cellAt(x-1, y)->type()==SOLID || cellAt(x, y)->type()==SOLID;
				if (borders_fluid_cell && !borders_solid_cell) {
						float left_pressure;
						float cell_pressure;
						if (cellAt(x-1, y)->type()==FLUID) {
							int left_id = cellAt(x-1, y)->id();
							left_pressure = p[left_id];
						} else if (cellAt(x-1, y)->type()==AIR) {
							left_pressure = 1;
						}
						if (cell->type()==FLUID) {
							int cell_id = cell->id();
							cell_pressure = p[cell_id];
						} else if (cell->type()==AIR) {
							cell_pressure = 1;
						}
						x_grad = cell_pressure - left_pressure;
				}			
			}
			// y gradient
			float y_grad = 0;
			if (cellAt(x, y-1) != NULL) {
				bool borders_fluid_cell = cellAt(x, y-1)->type()==FLUID || cellAt(x, y)->type()==FLUID;
				bool borders_solid_cell = cellAt(x, y-1)->type()==SOLID || cellAt(x, y)->type()==SOLID;

				if (borders_fluid_cell && !borders_solid_cell) {
						float below_pressure;
						float cell_pressure;
						if (cellAt(x, y-1)->type()==FLUID) {
							int below_id = cellAt(x, y-1)->id();
							below_pressure = p[below_id];

						} else if (cellAt(x, y-1)->type()==AIR) {
							below_pressure = 1;
						}
						if (cell->type()==FLUID) {
							int cell_id = cell->id();
							cell_pressure = p[cell_id];
						} else if (cell->type()==AIR) {
							cell_pressure = 1;
						}
						y_grad = cell_pressure - below_pressure;
				}
			}

			Eigen::Vector2d pressure_gradient(x_grad, y_grad);
			Eigen::Vector2d diff = t/(fluidDensity*this->cellSize())*pressure_gradient;
			cell->setU(cell->u()-diff);
		}
	}
}

double MacGrid::getDivergence(int x, int y)
{
	double u[2], ux, uy;
	u[0] = 0.0;
	u[1] = 0.0;
	ux = 0.0;
	uy = 0.0;

	int oX, oY;
	GridCell *cell, *n;

	oX = x;
	oY = y;
	cell = this->cellAt(oX, oY);
	if(cell != NULL && cell->type() != SOLID)
	{
		n = this->cellAt(oX - 1, oY);
		if(n != NULL && n->type() != SOLID)
			u[0] = cell->u()[0];
		n = this->cellAt(oX, oY - 1);
		if(n != NULL && n->type() != SOLID)
			u[1] = cell->u()[1];
	}

	oX = x + 1;
	oY = y;
	cell = this->cellAt(oX, oY);
	if(cell != NULL && cell->type() != SOLID)
		ux = cell->u()[0];

	oX = x;
	oY = y + 1;
	cell = this->cellAt(oX, oY);
	if(cell != NULL && cell->type() != SOLID)
		uy = cell->u()[1];

	return ux - u[0] + uy - u[1];
}

void MacGrid::extrapolateVelocity(int kcfl)
{
	//First set the layer of each fluid cell to 0 and -1 for all others
	GridCell *cell, *n0, *n1;
	GridCell *neighbors[4];
	for(int x = 0; x < this->_width_; ++x)
	{
		for(int y = 0; y < this->_height_; ++y)
		{

			cell = this->cellAt(x, y);
			if(cell->type() == FLUID)
				cell->setLayer(0);
			else
				cell->setLayer(-1);
		}
	}

	bool nInLayer;
	int lRange = max(2, kcfl);
	for(int layer = 1; layer <= lRange; ++layer)
	{
		for(int x = 0; x < this->_width_; ++x)
		{
			for(int y = 0; y < this->_height_; ++y)
			{
				cell = this->cellAt(x, y);
				if(cell->layer() == -1)
				{
					nInLayer = false;
					this->getNeighbors(x, y, neighbors);
					for(int i = 0; i < 4; ++i)
						if(neighbors[i] != NULL && neighbors[i]->layer() == layer - 1)
							nInLayer = true;

					if(nInLayer == true)
					{
						//set velocity components not bordering FLUID to the average of neighbors
						//in layer layer - 1
						n0 = this->cellAt(x - 1, y);
						if((n0 == NULL || n0->type() != FLUID) && cell->type() != FLUID)
						{
							n1 = this->cellAt(x + 1, y);
							if(n0 != NULL && n0->layer() == layer - 1)
							{
								if(n1 != NULL && n1->layer() == layer - 1)
								{
									cell->u()[0] = (n0->u()[0] + n1->u()[0]) / 2.0;
								}
								else
								{
									cell->u()[0] = n0->u()[0];
								}
							}
							else if(n1 != NULL && n1->layer() == layer - 1)
							{
								cell->u()[0] = n1->u()[0];
							}


						}

						n0 = this->cellAt(x, y - 1);
						if((n0 == NULL || n0->type() != FLUID) && cell->type() != FLUID)
						{
							n1 = this->cellAt(x, y + 1);
							if(n0 != NULL && n0->layer() == layer - 1)
							{
								if(n1 != NULL && n1->layer() == layer - 1)
								{
									cell->u()[1] = (n0->u()[1] + n1->u()[1]) / 2.0;
								}
								else
								{
									cell->u()[1] = n0->u()[1];
								}
							}
							else if(n1 != NULL && n1->layer() == layer - 1)
							{
								cell->u()[1] = n1->u()[1];
							}
						}

					}
				}
			}
		}
	}
}

void MacGrid::setSolidVelocities(void)
{
	GridCell *cell;

	for(int x = 0; x < this->_width_; ++x)
	{
		for(int y = 0; y < this->_height_; ++y)
		{
			cell = this->cellAt(x, y);

			if(cell->type() != FLUID && cell->type() != AIR)
				continue;

			if(x == 0 && cell->u()[0] < 0.0)
				cell->u()[0] = 0.0;
			if(y == 0 && cell->u()[1] < 0.0)
				cell->u()[1] = 0.0;
		}
	}
}

GridCell* MacGrid::cellAt(int i, int j) const
{
	if(i >= this->_width_ || j >= this->_height_)
		return NULL;
	else if(i < 0 || j < 0)
		return NULL;
	else
		return this->_cells_[i][j];
}

GridCell* MacGrid::cellAtWorldPos(double x, double y) const
{
	int i = floor(x / this->cellSize());
	int j = floor(y / this->cellSize());
	return cellAt(i, j);
}

void MacGrid::getNeighbors(int i, int j, GridCell** result) const
{
	result[0] = this->cellAt(i - 1, j);
	result[1] = this->cellAt(i, j - 1);
	result[2] = this->cellAt(i + 1, j);
	result[3] = this->cellAt(i, j + 1);
}

void MacGrid::serialize(const string path)
{
	ofstream f;
	f.open(path.c_str());

	GridCell *cell;
	Eigen::Vector2d pos, u;
	pos[0] = 0;
	pos[1] = 0;
	int index = 0;
	for(int y = 0; y < this->_height_; ++y)
	{
		pos[0] = 0;
		for(int x = 0; x < this->_width_; ++x)
		{
			cell = this->_cells_[x][y];
			u = cell->u();
			f << pos[0] << " " << pos[1] << " " << u[0] << " " << u[1] << " " <<
					this->cellSize() << " " << index++ << " " << cell->type() << " " <<
					cell->layer() << " " << cell->p() << "\n";
			pos[0] += this->cellSize();
		}
		pos[1] += this->cellSize();
	}

	f.close();
}

/*
 * private function definitions
 */

double MacGrid::getInterpolatedValue(double x, double y, int index) const
{
	int i = floor(x);
	int j = floor(y);

	double weights[4], vels[4];
	this->getInterpWeights(x, y, i, j, weights);
	this->getCellUComponents(i, j, index, vels);

	return weights[0] * vels[0] +
		   weights[1] * vels[1] +
		   weights[2] * vels[2] +
		   weights[3] * vels[3];
}

double MacGrid::getInterpolatedDiff(double x, double y, int index) const
{
	int i = floor(x);
	int j = floor(y);

	double weights[4], vels[4];
	this->getInterpWeights(x, y, i, j, weights);
	this->getCellDiffComponents(i, j, index, vels);
	// std::cout << "diffs: " << vels[0] << '\n';
	return weights[0] * vels[0] +
		   weights[1] * vels[1] +
		   weights[2] * vels[2] +
		   weights[3] * vels[3];
}

void MacGrid::getInterpWeights(double x, double y, int i, int j, double* result) const
{
    result[0] = (i + 1 - x) * (j + 1 - y);
    result[1] = (x - i) * (j + 1 - y);
    result[2] = (i + 1 - x) * (y - j);
    result[3] = (x - i) * (y - j);
}

double MacGrid::getCellU(int i, int j, int index) const
{
	GridCell *cell = this->cellAt(i, j);

	if(cell == NULL)
		return NON_EXISTENT_VEL;
	else
		return cell->u()[index];
}

double MacGrid::getCellUOld(int i, int j, int index) const
{
	GridCell *cell = this->cellAt(i, j);

	if(cell == NULL)
		return NON_EXISTENT_VEL;
	else
		return cell->oldU()[index];
}

void MacGrid::getCellUComponents(int i, int j, int index, double* result) const
{
	result[0] = this->getCellU(i, j, index);
	result[1] = this->getCellU(i + 1, j, index);
	result[2] = this->getCellU(i, j + 1, index);
	result[3] = this->getCellU(i + 1, j + 1, index);
}

void MacGrid::getCellDiffComponents(int i, int j, int index, double* result) const
{
	result[0] = this->getCellU(i, j, index) - this->getCellUOld(i, j, index);
	result[1] = this->getCellU(i + 1, j, index) - this->getCellUOld(i + 1, j, index);
	result[2] = this->getCellU(i, j + 1, index) - this->getCellUOld(i, j + 1, index);
	result[3] = this->getCellU(i + 1, j + 1, index) - this->getCellUOld(i + 1, j + 1, index);
}

int MacGrid::relabelFluidCells(void)
{
	int cellId = 0;
	GridCell* cell;
	for(int x = 0; x < this->_width_; ++x)
	{
		for(int y = 0; y < this->_height_; ++y)
		{
			cell = this->cellAt(x, y);
			if(cell->type() == FLUID)
				cell->setId(cellId++);
		}
	}

	return cellId;
}

void MacGrid::buildPressureMatrix(double t, double fluidDensity, double atmP)
{
	cout << "buildPressureMatrix: DONE!!!" << endl;
	int fluid_cell_count = this->relabelFluidCells();
	this->_b_->resize(fluid_cell_count);
	this->_A_->resize(fluid_cell_count, fluid_cell_count);
	this->_p_->resize(fluid_cell_count);
	float divergence;
	float air_cell_count;
	float non_solid_count;
	int cell_id;

	float h = this->_cellSize_;
	float density = fluidDensity;

	GridCell* neighbors[4];
	for (int x = 0; x < width(); ++x) {
		for (int y = 0; y < height(); ++y) {
			GridCell* cell = this->cellAt(x, y);

			if (cell->type() == FLUID) {
				cell_id = cell->id();
				air_cell_count = 0;
				non_solid_count = 0;
				this->getNeighbors(x, y, neighbors);

				for (int i = 0; i < 4; ++i) {
					if (neighbors[i] != NULL) {
						CellType type = neighbors[i]->type();
						if (type == AIR) {
							++air_cell_count;
							++non_solid_count;
						}
						if (type == FLUID) {
							++non_solid_count;
						}
					}
				}

				divergence = this->getDivergence(x, y);
				double b = density*h/t*divergence-air_cell_count*atmP;
				(*this->_b_)(cell_id) = b;
				this->_A_->insert(cell_id, cell_id) = -non_solid_count;
				
				for (GridCell* n : neighbors) {
					if (n != NULL) {
						CellType type = n->type();
						if (type == FLUID) {
							int index = n->id();
							this->_A_->insert(index, cell_id) = 1;
						}
					}
				}
			}
		}
	}
}
