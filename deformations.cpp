/**********************************************************************
* 
* deformations.cpp
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

#include <i3d/transform.h>
#include <i3d/FastLevelSet.h>
#include <i3d/morphology.h>
#include "deformations.h"
#include "perlin.h"

//----------------------------------------------------------------------------

void GenerateSpeedImage(i3d::Image3d<float> &speedImg, 
								double blob_size, 
								double smoothness, 
								float min_value, 
								float max_value, 
								float skewness)
{
	// generate Perlin noise 
	DoPerlin3D(speedImg, blob_size, smoothness, smoothness);
}

//----------------------------------------------------------------------------

void Deformation(i3d::Image3d<bool> &mask, 
					  const i3d::Image3d<bool> *domain, 
					  const i3d::Image3d<float> &speedImg, 
					  float deformation_level)
{
	 //+++++++++++
	 // PATCH
	 // If the deformed area is completely uniform (black/white), do nothing.
	 bool is_black(true), is_white(true);

	 if (domain == NULL)
	 {
		  for (size_t i=0; i<mask.GetImageSize(); i++)
		  {
				if (mask.GetVoxel(i))
					 is_black = false;
				else
					 is_white = false;
		  }
	 }
	 else
	 {
		  for (size_t i=0; i<mask.GetImageSize(); i++)
		  {
				if (domain->GetVoxel(i))
				{
					 if (mask.GetVoxel(i))
						  is_black = false;
					 else
						  is_white = false;
				}
		  }
	 }

	 if (is_black || is_white) // the mask is completely uniform
		  return; // leave this function without doing anything

	 // END OF PATCH
	 //+++++++++++

	 // mirror the speed image in order to have zero central differences on the border
	 i3d::Image3d<float> speedImgExt;
	 i3d::Mirror(speedImg, speedImgExt);


	 i3d::Vector3d<float> res = mask.GetResolution().GetRes();

	// we can use only one resolution component because the image is isotropic
	float inflation_impact = 1.0f;
	float curvature_impact = 0.5f;
	int radius = (int) floorf(0.4f * res.x); 

	try
	{
		i3d::fls::Routine<int> *lsSolver = i3d::fls::DefRoutine(speedImgExt, mask, domain, radius, inflation_impact, curvature_impact);

		lsSolver -> Initialize();

		int iteration = (int) floorf(lsSolver -> InterfaceSize() * deformation_level);

		lsSolver -> Iterate(iteration);
		lsSolver -> Mask(mask.begin(), i3d::fls::ResultCreator::InteriorAndInterface, true);

		delete lsSolver;

		// smooth the deformed object boundaries but protect the other objects 
		// from being modified by using mask image
		i3d::Image3d<bool> buffer;

		if (domain != NULL)
			i3d::LocalOpening(mask, buffer, *domain, i3d::nb3D_18);
		else
			i3d::Opening(mask, buffer, i3d::nb3D_18);

		mask = buffer;
	}
	catch (i3d::InternalException &e)
	{
		throw e.what.c_str();
	}
	catch (std::bad_alloc &)
	{
		throw "Deformation error: No memory."; 
	}
	catch (...)
	{
		throw "Deformation error: Unknown exception.";
	}
}

//----------------------------------------------------------------------------
