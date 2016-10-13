/*
 * GridCell.cpp
 *
 *  Created on: Jun 20, 2016
 *      Author: sflynn
 */

#include "GridCell.h"

GridCell::GridCell() : _type_(UNUSED), _id_(0), _layer_(0), _u_(0.0, 0.0), _tempU_(0.0, 0.0), _p_(0)
{

}

GridCell::~GridCell()
{

}

