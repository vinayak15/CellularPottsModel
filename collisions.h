/**********************************************************************
* 
* collisions.h
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
* Description: A collection of functions able to fit the newly generated
* cell into the cluster of already existing ones.
*
***********************************************************************/

#ifndef _COLLISIONS_H_
#define _COLLISIONS_H_

#include <i3d/image3d.h>

/**
 * This function finds all suitlable free position of an image mask
 * in the given space. 
 *
 * This function get two inputs: the mask of available area (0 - free,
 * 1 - occupied) and the mask of a new object (1 - foreground, 0 - background)
 * which we want to put into any free position of given available area.
 * These two images are submitted to the correlation function. All the possible
 * positions, where the object might be placed, correspond to zero correlation
 * coefficients.
 */
void LocateAllFreePositions(const i3d::Image3d<bool> &object,
									 const i3d::Image3d<bool> &ocuppiedSpace,
									 i3d::Image3d<bool> &freePositions);

/**
 * This function finds a suitable position of an image mask
 * in the given domain. It respects the clustering effect and
 * the gravity.
 *
 * This function get two inputs: the mask of available area (0 - free,
 * 1 - occupied) and the mask of a new object (1 - foreground, 0 - background)
 * which we want to put into any free position of given available area.
 * These two images are submitted to the correlation function. All the possible
 * positions, where the object might be placed, correspond to zero correlation
 * coefficients.
 * 
 * In order to guarantee one pixel gap between the object and the surrounding
 * objects the object is apriori submitted to one round of dilation.
 *
 * The positions are weighted with the distance from the already 
 * generated objects.
 *
 * This way, we generate the clustering effect of objects (cells).
 *
 * The detected position is the position of the center of the image which 
 * contains the observed object.
 *
 * return value: 
 * a) NULL ... position was not found (failure)
 * b) valid position 
 *
 * If a valid position is given, the memory location should be released manually
 * in the calling function.
 *
 */
i3d::Vector3d<float>* FindPosition(const i3d::Image3d<bool> &object,
											  const i3d::Image3d<bool> &space,
											  float cluster_effect = 25.0f,
											  bool influence_of_gravity = true);

/**
 * This function finds a suitable position of an (new) cell mask
 * in the given domain. It tries to find the position which is
 * the nearest to the positon of another (old) mask.
 *
 * This function get three inputs: the mask of available area (0 - free,
 * 1 - occupied), the mask of an old object (1 - foreground, 0 - background),
 * and the mask of a new object (1 - foreground, 0 - background)
 * We want to put a new object into any free position of given available area,
 * while it should be placed (if possible) to the same or near position of
 * the old one.
 * 
 * In order to guarantee one pixel gap between the object and the surrounding
 * objects the object is apriori submitted to one round of dilation.
 *
 * The positions are weighted with the distance from the already 
 * generated objects.
 *
 * This way, we generate the clustering effect of objects (cells).
 *
 * The detected position is the position of the center of the image which 
 * contains the observed object.
 *
 * return value: 
 * a) NULL ... position was not found (failure)
 * b) valid position 
 *
 * If a valid position is given, the memory location should be released manually
 * in the calling function.
 *
 */
i3d::Vector3d<float> *FindNewPosition(const i3d::Image3d<bool> &oldMask,
												  const i3d::Image3d<bool> &newMask,
												  const i3d::Image3d<bool> &space);

#endif
