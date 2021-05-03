
/**********************************************************************
*
* boundary.h
*
* This file is part of VesselGen(3D)
*
* Copyright (C) 2016 -- Centre for Biomedical Image Analysis (CBIA)
* http://cbia.fi.muni.cz/
*
* VesselGen is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* VesselGen is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with VesselGen. If not, see <http://www.gnu.org/licenses/>.
*
* Author: David Svoboda
*
* Description: Functions providing the manipulation with coordinates
* located near the image boundaries.
*
***********************************************************************/

#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_


/**
 * Apply periodic/non-periodic boundary condition check
 * on voxel coordinate. If not valid, the function return false.
 * If valid it returns true. In case of the periodic bounaries are
 * valid, the coordiates are change to be located in the first-period
 * within the domain.
 *
 * @param[in,out] coor		inspected coordinate
 * @param[in] sz	size 		of the image
 */
bool ValidateCoords(i3d::Vector3d<int> &coor, i3d::Vector3d<size_t> sz);


#endif
