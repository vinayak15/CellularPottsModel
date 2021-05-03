/**********************************************************************
* 
* perlin.h
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
* Description: An I3D wrapper for standard Perlin noise functions.
*
***********************************************************************/

#ifndef _PERLIN_H_
#define _PERLIN_H_

#include <i3d/image3d.h>

/**
 * Generate Perlin noise into the given image.
 */
void DoPerlin3D(i3d::Image3d<float> &fimg, 
					 double var, 
					 double alpha = 8, 
					 double beta = 4, 
					 int n = 6);


#endif
