/**********************************************************************
*
* edge_detection.h
*
* Author: Peter Kováč, David Svoboda
*
* Description: Functions to handle detecting and updating collection of 
* "edge voxels".
*
***********************************************************************/

#ifndef EDGE_DETECTION_H
#define EDGE_DETECTION_H

#include <i3d/image3d.h>
#include <i3d/neighbours.h>
//#include "cpm.h"
#include "settings.h"
#include <limits>


/**
  * Find "edge voxels" of image
  * Edge voxel have one of two types of connections with 
  * some voxel from his neighbourhood:
  *			A) MEDIUM - CELL 
  *			B) 2 DIFFERENT CELLS
  * 
  * Function is called only at initialization
  * 
  * @param[in]	scene		3D lattice with cell population
  * @param[in]	nbh			neighbourhood of voxel
  * @param[out] edgeSet		set with coordinates of edge voxels
 */
/*void EdgeDetection(const i3d::Image3d<i3d::GRAY16> &scene,
						 const i3d::Neighbourhood &nbh,
						 std::set<i3d::Vector3d<int> > &edgeSet );*/

/**
  * Decide if voxels in neighbourhood of "updated voxel"
  *	still are "edge voxels" and update edgeSet accordingly
  * 
  * @param[in]	scene			3D lattice with cell population
  * @param[in]	nbh			neighbourhood of voxel
  * @param[out] edgeSet		set with coordinates of edge voxels
  * @param[in]	upd_coor		updated voxel
 */
/*void UpdateEdge(const i3d::Image3d<i3d::GRAY16> &scene,
					 const i3d::Neighbourhood &nbh,
					 std::set<i3d::Vector3d<int> > &edgeSet,
					 i3d::Vector3d<int> &upd_voxel,
					 bool periodic);*/

const size_t nonEdge = std::numeric_limits<size_t>::max();

class Edges
{
  private: 
			 size_t *maskOfEdges, *shakedArrayOfEdges;
			 size_t numOfEdges;

  public: Edges(const i3d::Image3d<i3d::GRAY16> &scene,
					 const i3d::Neighbourhood &nbh);
			 ~Edges();

			 void Update(const i3d::Image3d<i3d::GRAY16> &scene,
							 const i3d::Neighbourhood &nbh, 
							 i3d::Vector3d<int> &upd_voxel ,const float *);

			 size_t size() const { return numOfEdges; };
			 size_t GetIndex(size_t i) const { return shakedArrayOfEdges[i]; };
};

#endif	/* EDGE_DETECTION_H */

