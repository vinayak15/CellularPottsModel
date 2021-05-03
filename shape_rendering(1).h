/**********************************************************************
* 
* shape_rendering.h
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
* Description: Functions for generating the basic rough cell shape plus
* inserting it into another image containing the already generated cell
* population.
*
***********************************************************************/

#ifndef SHAPE_RENDERING_H
#define SHAPE_RENDERING_H

#include <i3d/image3d.h>

#include "settings.h"

/**
 * Generate a black cell mask with a given resolution.
 * The image variable 'mask' needn't be initialized as
 * the function will take care of it.
 *
 * \param[in] res	required resolution of newly created image
 * \param[out]	mask	newly created image with a cell mask inside
 * \param[in] cellDiameter mean diameter (in microns) of the cell mask that should be generated
 */
void RenderBlankCellMask(const i3d::Vector3d<float> &res,
								 i3d::Image3d<bool> &mask,
								 float cellDiameter,
								 float shape_variance = 16.0);

/**
 * Take an image a draw an generic ellipsoid inside.
 *
 * \param[in,out] img	prepared and allocated free image buffer
 * \param[in] center	position of center of ellipsoid
 * \param[in] radius size of semiaxes of ellipsoid
 * \param[in] angle_xy angle of rotation of ellipsoid in xy plane 
 * \param[in] angle_xz angle of rotation of ellipsoid in xz plane
 * \param[in] foregoundIntensity intensity/color of ellipsoid voxels
 */
template <class VOXEL> void RenderEllipsoid(
		i3d::Image3d<VOXEL> &img,
		const i3d::Vector3d<float> &center,
		const i3d::Vector3d<float> &radius,
		const float angle_xy,
		const float angle_xz,
		const VOXEL foregroundIntensity);

/** 
 * Merge the binary mask with the labeled image
 *
 * \param[in,out] img labeled image
 * \param[in] mask merged binary image
 */
template <class T> bool AddNewLabel(i3d::Image3d<T> &img, 
												const i3d::Image3d<bool>& mask);

/**
 * Copy out the mask from the labeled image + remember offset
 *
 * \param[in] img original image with labeled objects
 * \param[in] label selected label
 * \param[out] mask object with given label extracted from the original image
 */
template <class T> void GetLabelMask(const i3d::Image3d<T> &img,
												 T label,
												 i3d::Image3d<bool> &mask);
#endif
