/**********************************************************************
* 
* deformations.h
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
* Author: Martin Maska, David Svoboda
* 
* Description: The topology-preserving deformation of binary objects 
* using a fast level set-like algorithm.
*
***********************************************************************/

#ifndef _DEFORMATION_H_
#define _DEFORMATION_H_

#include <i3d/image3d.h>

/** The topology-preserving deformation of binary objects using a fast level set-like algorithm. The boundary of objects 
  * in the mask image is taken as the initial interface that is deformed under a force field given by the speed image. 
  * The last parameter determines the level of deformation, i.e. how much the initial interface is deformed. The second
  * parameter defines the computational domain, i.e. the area in which the interface can be evolved. It is recommended
  * to set this parameter to NULL if a full image deformation is required. Deformed objects are saved back to the mask 
  * image. Note that if a domain is specified then only those points belonging to the domain are saved back, the remaining
  * ones are left without any change. */
void Deformation(i3d::Image3d<bool> &mask, 
					  const i3d::Image3d<bool> *domain, 
					  const i3d::Image3d<float> &speedImg, 
					  float deformation_level);

/** Generate a speed image. */
void GenerateSpeedImage(i3d::Image3d<float> &speedImg, 
								double blob_size, 
								double smoothness = 8, 
								float min_value = 0.0f, 
								float max_value = 1.0f, 
								float skewness = 1.0f);

#endif // _DEFORMATION_
