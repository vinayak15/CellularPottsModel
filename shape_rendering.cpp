/**********************************************************************
* 
* shape_rendering.cpp
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


#include <math.h>
#include <i3d/FastLevelSet.h>
#include <i3d/image3d.h>
#include <i3d/descriptors/shape.h>
#include <i3d/morphology.h>

#include "shape_rendering.h"
#include "rnd_generators.h"
#include "deformations.h"

//----------------------------------------------------------------
void RenderBlankCellMask(const i3d::Vector3d<float> &res,
								 i3d::Image3d<bool> &mask,
								 float cellDiameter,
								 float shape_variance)
{
	 // prepare the initial estimate of object dimensions
	 i3d::Vector3d<float> r;
	 r.x = GetRandomGauss (cellDiameter/2.0f, (shape_variance/100.0f) * cellDiameter/2.0f);
	 r.y = cellDiameter/2.0f + (cellDiameter/2.0f - r.x);
	 r.z = GetRandomGauss (cellDiameter/3.0f, (shape_variance/100.0f) * cellDiameter/3.0f);

	 float xy_rotation = GetRandomGauss(0.0f, 45.0f);
	 float xz_rotation = GetRandomGauss(0.0f, 10.0f);

	 i3d::Vector3d<size_t> sz;
	 float max_r = std::max(r.x, std::max(r.y, r.z));
	 i3d::Vector3d<float> centre(max_r + 1.0f), dims = centre * 2.0f; 
	 
	 sz = i3d::MicronsToPixels(dims, res);
	 mask.MakeRoom(sz); // plus some boundary
	 mask.SetResolution(res);
	 mask.Clear();

	 DEBUG_REPORT("allocating size: " << sz);
	 DEBUG_REPORT("with resolution: " << res);

	 // create the cellMask object
	 RenderEllipsoid (mask, centre, r, xy_rotation, xz_rotation, true);
	 DEBUG_REPORT("ellipsoid semiaxes: " << r);

	 // slightly deform the initial object
	 i3d::Image3d<float> speedImg;
	 speedImg.CopyMetaData(mask);
	
	 GenerateSpeedImage(speedImg, (r.x + r.y + r.z)/4.5f);
	 float deformationLevel = 0.7f;
	 Deformation(mask, NULL, speedImg, deformationLevel); 
}

//----------------------------------------------------------------

template <class T> bool AddNewLabel(i3d::Image3d<T> &img, 
												const i3d::Image3d<bool>& mask)
{
	 if (img.GetResolution().GetRes() != mask.GetResolution().GetRes())
		  throw "Resolutions do not match.";

	 if (mask.GetMaxValue() == false) // nothing to add
	 {
		  DEBUG_REPORT("exiting -- nothing to add");
		  return false;
	 }

	 i3d::Offset diff_offset = mask.GetOffset() - img.GetOffset();
	 i3d::Vector3d<size_t> pos = i3d::MicronsToPixels(diff_offset, img.GetResolution());

	 size_t x,y,z;

	 T new_label = img.GetMaxValue() + 1;

	  for (size_t i=0; i<mask.GetImageSize(); i++)
	  {
			if (mask.GetVoxel(i))
			{
				 x = pos.x + mask.GetX(i);
				 y = pos.y + mask.GetY(i);
				 z = pos.z + mask.GetZ(i); 
					 
				 img.SetVoxel(x,y,z,new_label);
			}
	  }

	 return true;
}

//----------------------------------------------------------------

template <class T> void GetLabelMask(const i3d::Image3d<T> &img,
												 T label,
												 i3d::Image3d<bool> &mask)
{
	 size_t x,y,z;
	 size_t min_x, min_y, min_z, max_x, max_y, max_z;

	 // initialize the variables
	 min_x = img.GetSizeX();
	 min_y = img.GetSizeY();
	 min_z = img.GetSizeZ();

	 max_x = max_y = max_z = 0;

	 // find the voxels with given label
	 for (size_t i=0; i<img.GetImageSize(); i++)
	 {
		  if (img.GetVoxel(i) == label)
		  {
				x = img.GetX(i);
				y = img.GetY(i);
				z = img.GetZ(i); 
				
				min_x = (min_x > x) ? x : min_x;
				min_y = (min_y > y) ? y : min_y;
				min_z = (min_z > z) ? z : min_z; 
				
				max_x = (max_x < x) ? x : max_x; 
				max_y = (max_y < y) ? y : max_y;
				max_z = (max_z < z) ? z : max_z;
		  }
	 }

	 //enlarge the mask so that the object, HOPEFULLY, does not touch the border
	 //this is required for proper work of ResetBoundaryPoints() function
	 if (min_x > 2) min_x-=2; else min_x=0;
	 if (min_y > 2) min_y-=2; else min_y=0;
	 if (min_z > 2) min_z-=2; else min_z=0;
	 if (max_x < img.GetSizeX()-2) max_x+=2; else max_x=img.GetSizeX()-1;
	 if (max_y < img.GetSizeY()-2) max_y+=2; else max_y=img.GetSizeY()-1;
	 if (max_z < img.GetSizeZ()-2) max_z+=2; else max_z=img.GetSizeZ()-1;

	 // form the mask boundaries
	 i3d::Vector3d<size_t> label_offset(min_x, min_y, min_z);
	 i3d::Vector3d<size_t> label_size(max_x - min_x + 1,
										  max_y - min_y + 1,
										  max_z - min_z + 1);

	 mask.SetResolution(img.GetResolution());
	 mask.SetOffset(img.GetOffset() + 
						 PixelsToMicrons(label_offset, img.GetResolution()));
	 mask.MakeRoom(label_size);
	 mask.Clear();

	 for (x=0; x<label_size.x; x++)
		  for (y=0; y<label_size.y; y++)
				for (z=0; z<label_size.z; z++)
				{
					 if (img.GetVoxel(label_offset.x + x,
											label_offset.y + y,
											label_offset.z + z) == label)
					 {
						  mask.SetVoxel(x,y,z,true);
					 }
				}
}

//----------------------------------------------------------------


inline bool IsInEllipsoid(float x, float y, float z,
		float sx, float sy, float sz,
		float rx, float ry, float rz,
		float rotate_xy, float rotate_xz)
{
	 // first proceed fast rough check
	if (SQR(x-sx) + SQR(y-sy) + SQR(z-sz) > 
		 SQR(std::max(std::max(rx, ry), rz)))
		  return false;

	float x1, x2, x3, y1, y2, y3, z1, z2, z3;

	// shift to center
	x1 = x - sx;
	y1 = y - sy;
	z1 = z - sz;

	// xy-plane rotation
	x2 = float(x1 * cos(DEG_TO_RAD(rotate_xy)) 
				  - y1 * sin(DEG_TO_RAD(rotate_xy)));

	y2 = float(x1 * sin(DEG_TO_RAD(rotate_xy)) 
				  + y1 * cos(DEG_TO_RAD(rotate_xy)));

	z2 = z1;

	// xz-plane rotation
	x3 = float(x2 * cos(DEG_TO_RAD(rotate_xz)) 
				  - z2 * sin(DEG_TO_RAD(rotate_xz)));

	y3 = y2;

	z3 = float(x2 * sin(DEG_TO_RAD(rotate_xz)) 
				  + z2 * cos(DEG_TO_RAD(rotate_xz)));

	return ((SQR(x3/rx) + SQR(y3/ry) + SQR(z3/rz)) <= 1.0f);
}


//----------------------------------------------------------------

template <class VOXEL> void RenderEllipsoid(
		i3d::Image3d<VOXEL> &img,
		const i3d::Vector3d<float> &center,
		const i3d::Vector3d<float> &radius,
		const float angle_xy,
		const float angle_xz,
		const VOXEL foregroundIntensity)
{
	 i3d::Vector3d<float> res = img.GetResolution().GetRes();
	
	for (size_t i=0; i<img.GetImageSize(); i++)
	{
		if (IsInEllipsoid(
					img.GetX(i)/res.x, img.GetY(i)/res.y, img.GetZ(i)/res.z,
					center.x, center.y, center.z,
					radius.x, radius.y, radius.z,
					angle_xy, angle_xz))
		{
			img.SetVoxel(i, foregroundIntensity);
		}
	}
}

/*************************************************************************/
/* explicit instantiations */
template 
void RenderEllipsoid(i3d::Image3d<bool> &,
							const i3d::Vector3d<float> &,
							const i3d::Vector3d<float> &,
							const float angle_xy,
							const float angle_xz,
							const bool foreground_intensity);

template 
void RenderEllipsoid(i3d::Image3d<i3d::GRAY8> &,
							const i3d::Vector3d<float> &,
							const i3d::Vector3d<float> &,
							const float angle_xy,
							const float angle_xz,
							const i3d::GRAY8 foreground_intensity);

template 
void RenderEllipsoid(i3d::Image3d<i3d::GRAY16> &,
							const i3d::Vector3d<float> &,
							const i3d::Vector3d<float> &,
							const float angle_xy,
							const float angle_xz,
							const i3d::GRAY16 foreground_intensity);

template 
void RenderEllipsoid(i3d::Image3d<size_t> &,
							const i3d::Vector3d<float> &,
							const i3d::Vector3d<float> &,
							const float angle_xy,
							const float angle_xz,
							const size_t foreground_intensity);


template 
void GetLabelMask(const i3d::Image3d<i3d::GRAY8> &img,
						i3d::GRAY8 label,
						i3d::Image3d<bool> &mask);
template 
void GetLabelMask(const i3d::Image3d<i3d::GRAY16> &img,
						i3d::GRAY16 label,
						i3d::Image3d<bool> &mask);

template bool AddNewLabel(i3d::Image3d<i3d::GRAY16> &img, 
								  const i3d::Image3d<bool>& mask);
