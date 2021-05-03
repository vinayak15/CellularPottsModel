/**********************************************************************
* 
* collisions.cc
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

#include <float.h>
#include <i3d/DistanceTransform.h>
#include <i3d/morphology.h>
#include <i3d/se.h>
#include <i3d/convolution.h>
#include <i3d/descriptors/shape.h>

#include "collisions.h"
#include "settings.h"

#define SMALL_VALUE 0.5f
#define TOLERANCE_IN_MICRONS 0.4f

//---------------------------------------------------------------------------
template <class T>
bool GetLastNonZeroPosition(const i3d::Image3d<T> &img, i3d::Vector3d<size_t> &v)
{
	 for (size_t i = img.GetImageSize()-1; i>=0; i--)
	 {
		  if (img.GetVoxel(i) != T(0))
		  {
				v = img.GetPos(i);
				return true;
		  }
	 }

	 return false;
}

//---------------------------------------------------------------------------

void LocateAllFreePositions(const i3d::Image3d<bool> &object,
									 const i3d::Image3d<bool> &ocuppiedSpace,
									 i3d::Image3d<bool> &freePositions)
{
	 i3d::Image3d<bool> objectCopy(object), spaceDilated;
	 i3d::Image3d<float> correlationResults;

	 DEBUG_REPORT("dilation of space image");

	 // dilate the space to guarantee there will be a gap between
	 // the generated objects
	 i3d::Dilation(ocuppiedSpace, spaceDilated, i3d::nb3D_26);

	 DEBUG_REPORT("drawing border line");

	 // draw a border line in the space to prevent the possible object
	 // positions from being located at the border. 
	 for (size_t i=0; i<spaceDilated.GetImageSize(); i++)
	 {
		  i3d::Vector3d<size_t> coordinates = spaceDilated.GetPos(i);

		  if (spaceDilated.OnBorder(coordinates))
		  {
				spaceDilated.SetVoxel(i, true);
		  }
	 }

	 DEBUG_REPORT("performing convolution");

	 // apply convolution; as we flip the image in advance we apply
	 // the correlation instead of convolution
	 objectCopy.Flip(true, true, true);
	 i3d::Convolution<double>(spaceDilated, 
									  objectCopy, 
									  correlationResults, 
									  false);

	 // Find the weights by distance transform:
	 freePositions.CopyMetaData(correlationResults);

	 DEBUG_REPORT("thresholding convolution result");

	 // Mask the results of distance transform with the correlation results
	 for (size_t i=0; i<correlationResults.GetImageSize(); i++)
	 {
		  // remove the positions, where the correlation results
		  // give some non-zero value
		  if (correlationResults.GetVoxel(i) < SMALL_VALUE)
				freePositions.SetVoxel(i, true);
		  else
				freePositions.SetVoxel(i, false);
	 }
}

//---------------------------------------------------------------------------
i3d::Vector3d<float> *FindNewPosition(const i3d::Image3d<bool> &oldMask,
												  const i3d::Image3d<bool> &newMask,
												  const i3d::Image3d<bool> &space)
{
	 i3d::Image3d<bool> free_space;
	 free_space.CopyMetaData(space);

	 LocateAllFreePositions(newMask, space, free_space);

	 // localize centre of mass for old mask (absolute position in space)
	 i3d::Shape<bool> oldShape;
	 oldShape.SetImage(&oldMask);
	 oldShape.SetLabel(true);

	 i3d::Vector3d<float> oldMaskCentreAbs = 
				i3d::PixelsToMicrons(oldShape.Centroid(),
											oldMask.GetResolution());
	 oldMaskCentreAbs += oldMask.GetOffset();

	 // localize centre of mass for new mask (absolute position within the subimage)
	 i3d::Shape<bool> newShape;
	 newShape.SetImage(&newMask);
	 newShape.SetLabel(true);

	 i3d::Vector3d<float> newMaskCentreRel = 
				i3d::PixelsToMicrons(newShape.Centroid(),
											newMask.GetResolution());

	 // compute the half size of new image in microns
	 i3d::Vector3d<float> newImgCentre = 
				i3d::PixelsToMicrons(newMask.GetSize(),
											newMask.GetResolution()) / 2.0f;

	 // find the position where the distance between the old mask and the new mask  
	 // is the shortest
	 float minDistance = FLT_MAX;
	 size_t index = 0;

	 for (size_t i=0; i<free_space.GetImageSize(); i++)
	 {
		  // check only the position with non-zeros values, i.e.
		  // the places suitable for nesting the new object
		  if (free_space.GetVoxel(i) == true)
		  {
				// fpos = absolute location of centre of image within the global space
				i3d::Vector3d<float> fpos = 
						  i3d::PixelsToMicrons(free_space.GetPos(i),
													  free_space.GetResolution()) + 
						  free_space.GetOffset();

				// newMaskCentreAbs = absolute location of a new mask within the global
				// space
				i3d::Vector3d<float> newMaskCentreAbs =
						  fpos + newMaskCentreRel - newImgCentre;

				// distance between centre of mass of old cell mask and the newly 
				// created cell mask
				float currDistance = i3d::Norm2(newMaskCentreAbs - oldMaskCentreAbs);

				if (currDistance < minDistance)
				{
					 index = i;
					 minDistance = currDistance;
				}
		  }
	 }

	 if (minDistance == FLT_MAX)
		  return NULL;

	 // do not forget to free this variable in the calling function!
	 i3d::Vector3d<float> *position = new i3d::Vector3d<float>;

	 *position = i3d::PixelsToMicrons(free_space.GetPos(index), 
												 free_space.GetResolution());

	 *position += free_space.GetOffset();

	 return position;
}


//---------------------------------------------------------------------------

i3d::Vector3d<float> *FindPosition(const i3d::Image3d<bool> &object, 
											  const i3d::Image3d<bool> &space,
											  float cluster_effect,
											  bool influence_of_gravity)
{
	 i3d::Image3d<bool> free_space;
	 free_space.CopyMetaData(space);

	 LocateAllFreePositions(object, space, free_space);

	 // if required impose the gravity to the objects
	 if (influence_of_gravity)
	 {
		  i3d::Vector3d<float> res = free_space.GetResolution().GetRes();
		  size_t tolerance_in_pixels = 
					 (size_t) floorf(TOLERANCE_IN_MICRONS * res.z);

		  i3d::Vector3d<size_t> pos;

		  DEBUG_REPORT("tolerance in pixels: " << tolerance_in_pixels);

		  if (!GetLastNonZeroPosition(free_space, pos))
		  {
				throw ERROR_REPORT("No more space for new object");
		  }

		  DEBUG_REPORT("last non zero position: " << pos);
		  size_t spread = 
			 (pos.z >= tolerance_in_pixels) ? (pos.z - tolerance_in_pixels) : 0;

		  for (size_t i=0; i<free_space.GetImageSize(); i++)
		  {
				if (free_space.GetZ(i) < spread)
					 free_space.SetVoxel(i, false);
		  }
	 }

	 // Respect the cluster effect when generating new positions.
	 bool cluster = (((float)rand())/RAND_MAX) < (cluster_effect/100.0f);

	 // Apply cluster effect only if there exists at least one object. Otherwise,
	 // it has no sense.
	 if (cluster && (space.GetMaxValue() == true))
	 {
		  DEBUG_REPORT("cluster effect: yes");
		  i3d::Image3d<float> distance_map;

		  distance_map.CopyMetaData(space);

		  for (size_t i=0; i<distance_map.GetImageSize(); i++)
		  {
				if (space.GetVoxel(i) == true)
					 distance_map.SetVoxel(i, 1.0f);
				else
					 distance_map.SetVoxel(i, 0.0f);
		  }

		  i3d::EDM(distance_map, 0, 0.0f, false);

		  for (size_t i=0; i<distance_map.GetImageSize(); i++)
		  {
				if (free_space.GetVoxel(i) == false)
					 distance_map.SetVoxel(i, FLT_MAX);
		  }

		  float min_distance = distance_map.GetMinValue();

		  for (size_t i=0; i<distance_map.GetImageSize(); i++)
		  {
				if (free_space.GetVoxel(i) == true)
				{
					 float dist_value = distance_map.GetVoxel(i);

					 if ((dist_value < min_distance) ||
						  (dist_value > min_distance + TOLERANCE_IN_MICRONS))
					 {
						  free_space.SetVoxel(i, false);
					 }
				}
		  }
	 }

	 // check the number of free positions
	 size_t num_of_positions = 0;

	 for (size_t i=0; i<free_space.GetImageSize(); i++)
	 {
		  if (free_space.GetVoxel(i) == true)
				num_of_positions++;
	 }

	DEBUG_REPORT("sum of free positions: " << num_of_positions);

	 // If no suitable position was found, leave the function with no success.
	 if (num_of_positions == 0)
	 {
		  return NULL;
	 }

	 // Find the suitable position 
	 size_t barrier = (size_t) ((((float)rand()) / RAND_MAX)*num_of_positions);
	 size_t index = 0, shift = 0;

	 while (shift <= barrier)
	 {
		  shift += (free_space.GetVoxel(index)) ? 1 : 0;
		  index++;
	 }

	 // do not forget to free this variable in the calling function!
	 i3d::Vector3d<float> *position = new i3d::Vector3d<float>;

	 *position = i3d::PixelsToMicrons(free_space.GetPos(index), 
												 free_space.GetResolution());

	 *position += free_space.GetOffset();

	 DEBUG_REPORT("searching finished");

	 return position;
}

