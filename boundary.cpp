/**********************************************************************
*
* boundary.cpp
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

* Authors: David Svoboda,
*          Peter Kováč
*
* Description: Functions providing the manipulation with coordinates
* located near the image boundaries.
*
***********************************************************************/

#include <i3d/image3d.h>
#include "boundary.h"

bool ValidateCoords(i3d::Vector3d<int> &coor, i3d::Vector3d<size_t> sz)
{
	 // z-coordiate is always limited
	 if ((coor.z < 0) || (coor.z >= (int) sz.z))
		  return false; 

	 coor.x = (coor.x + (int) sz.x) % ((int) sz.x);
	 coor.y = (coor.y + (int) sz.y) % ((int) sz.y);


#ifdef BOUNDARY_LIMITED
	 if ((coor.x < 0) || (coor.x >= (int) sz.x))
				return false;

    if ((coor.y < 0) || (coor.y >= (int) sz.y))
				return false;
#endif

	 return true;
}

//-----------------------------------------------------------------------

