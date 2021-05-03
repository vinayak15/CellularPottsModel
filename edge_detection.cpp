#include "edge_detection.h"
#include "settings.h"
#include "boundary.h"

#include <i3d/neighbours.h>

//-----------------------------------------------------------------------

Edges::Edges(const i3d::Image3d<i3d::GRAY16> &scene,
				 const i3d::Neighbourhood &nbh)
{
	numOfEdges = 0;
	shakedArrayOfEdges = new size_t[scene.GetImageSize()];
	maskOfEdges = new size_t[scene.GetImageSize()];

	for (size_t i=0; i<scene.GetImageSize(); i++)
	{
		 maskOfEdges[i] = nonEdge;
	}

	DEBUG_REPORT("Searching for all edge voxels in the image ...");
	std::vector<const i3d::GRAY16 *> winNbh;
	
	// For all voxels (except borders) check if voxel has in 
	// neighbourhood voxel with another value(it means that voxel creates edge)
	// if test will pass, voxel's coordinates are added into edge_set
	size_t i,j;
	for (i=0, j=0; i < scene.GetImageSize(); i++)
	{
		if (scene.GetVoxel(i) != ID_ECM)
		{
			if (j==0)
			{
				// get neighbourhood of voxel
				j = i3d::GetWindow(scene, 
								scene.GetX(i), 
								scene.GetY(i), 
								scene.GetZ(i), 
								nbh, winNbh); 
			}
			else
			{
				 i3d::MoveWindow(winNbh);
				 j--;
			}
		
			bool edgeFound = false;

			for (unsigned int k = 1; (k < winNbh.size()) && (!edgeFound); k++) 
			{
				 if ((*winNbh[k] != ID_ECM) && (scene.GetVoxel(i) != *winNbh[k])) 
				 {
					  // If current voxel's value != neighbour voxel's value,
					  // neighbour voxel's value must be DIFFERENT from ECM	
					  // mark index of current voxel as an edge. The mask image
					  // points to array of edges and vice versa
					  maskOfEdges[i] = numOfEdges;
					  shakedArrayOfEdges[numOfEdges] = i;
					  numOfEdges++;

					  edgeFound  = true;
				 }
			}
		}
		else // just skip the ECM values 
		{
			 if (j>0)
			 {
				  i3d::MoveWindow(winNbh);
				  j--;
			 }
		}
	}
	DEBUG_REPORT("All edge voxels were detected.");
	DEBUG_REPORT("Number of edge voxels: " << numOfEdges);


	i3d::Image3d<i3d::GRAY16> imgEdges;
	imgEdges.CopyMetaData(scene);

	for (size_t i=0; i<imgEdges.GetImageSize(); i++)
	{
		 if (maskOfEdges[i] == nonEdge)
			  imgEdges.SetVoxel(i, 0);
		 else
			  imgEdges.SetVoxel(i, 65000);
	}
	imgEdges.SaveImage("_hrany.ics");

}

//-----------------------------------------------------------------------
Edges::~Edges()
{
	 delete [] shakedArrayOfEdges;
	 delete [] maskOfEdges;
}
//-----------------------------------------------------------------------
/*
void EdgeDetection(const i3d::Image3d<i3d::GRAY16> &scene,
						 const i3d::Neighbourhood &nbh,
						 std::set<i3d::Vector3d<int> > &edgeSet)
{
	DEBUG_REPORT("Finding all edge voxels in image...");
	std::vector<const i3d::GRAY16 *> winNbh;
	
	// for all voxels (except borders) check if voxel has in 
	// neighbourhood voxel with another value(it means that voxel creates edge)
	// if test will pass, voxel's coordinates are added into edge_set
	for (size_t i = 0; i < scene.GetImageSize(); i++)
	{
		if(scene.GetVoxel(i) == ID_ECM)
		{
			continue;
		}
		//get neighbourhood of voxel
		i3d::GetWindow(scene, scene.GetX(i), scene.GetY(i),
								scene.GetZ(i), nbh, winNbh);
		for (unsigned int j = 1; j < winNbh.size(); j++) {
			if (*winNbh[j] != ID_ECM && 
				scene.GetVoxel(i) != *winNbh[j]) {
				//if current voxel's value != neighbour voxel's value,
				//neighbour voxel's value must be DIFFERENT from ECM	
				//add coordinates of current voxel to edge_set
				edgeSet.insert(scene.GetPos(i));
				break;
			}
		}		
	}
	DEBUG_REPORT("All edge voxels were found.");
}
*/
//-----------------------------------------------------------------------

void Edges::Update(const i3d::Image3d<i3d::GRAY16> &scene,
	const i3d::Neighbourhood &nbh,
	i3d::Vector3d<int> &updVoxel, const float * imgCellIDs)
{
	// For storing values of voxel's neighbourhood
	std::vector<const float *> winNbh;

	// Iterate through all neighbours of updVoxel (including this voxel)
	for (size_t i = 0; i < nbh.size(); i++)
	{
		i3d::Vector3d<int> currVoxel(updVoxel + nbh.offset[i]);
		size_t currIndex = scene.GetIndex(currVoxel);

		// Handle boundary condition
		if (!ValidateCoords(currVoxel, scene.GetSize()))
		{
			continue;
		}

		if (imgCellIDs[currIndex] == ID_ECM)
		{
			continue;
		}

		// Get the values of currVoxel's neighbourhood
		i3d::Neighbourhood rnb;
		size_t tmp = GetNbh(scene, currVoxel.x, currVoxel.y, currVoxel.z, nbh, rnb);
		winNbh.resize(nbh.size());
		i3d::VectContainer::const_iterator off;
		int q = 0;
		for (off = nbh.offset.begin(); off != nbh.offset.end(); ++off)
		{
			winNbh[q++] = &imgCellIDs[currVoxel.x + off->x + (currVoxel.y + off->y)*scene.GetSizeX() + (currVoxel.z + off->z)*scene.GetSizeX() * scene.GetSizeY()];
		}
		bool isEdgeVoxel = false;

		for (size_t j = 1; (j < winNbh.size()) && (!isEdgeVoxel); j++)
		{
			// analyzing if currVoxel is "EDGE VOXEL"
			// compare with neighbourhood values
			// (nbh value NOT ECM and nbh, currVoxel value are DIFFERENT) ->
			//	-> currVoxel is "EDGE VOXEL"
			if ((*winNbh[j] != ID_ECM) &&
				(imgCellIDs[currIndex] != *winNbh[j]))
			{
				// If it was not an edge before, mark it as an edge voxel.
				if (maskOfEdges[currIndex] == nonEdge)
				{
					maskOfEdges[currIndex] = numOfEdges;
					shakedArrayOfEdges[numOfEdges] = currIndex;
					numOfEdges++;
				}

				isEdgeVoxel = true;
			}
		}

		// If voxel doesn't form an edge, unmark it if it was and edge voxel
		// befroe the change.
		if (!isEdgeVoxel)
		{
			if (maskOfEdges[currIndex] != nonEdge) // Formerly, it was an edge.
												   // We need to properly remove it.
			{
				if (currIndex == (numOfEdges - 1)) // It is the last in the array
				{
					maskOfEdges[currIndex] = nonEdge;
					numOfEdges--;
				}
				else // It appears to be located somewhere in the array of edges
				{
					maskOfEdges[shakedArrayOfEdges[numOfEdges - 1]] = maskOfEdges[currIndex];
					shakedArrayOfEdges[maskOfEdges[currIndex]] = shakedArrayOfEdges[numOfEdges - 1];
					maskOfEdges[currIndex] = nonEdge;
					numOfEdges--;
				}
			}
		}
	}
}


//-----------------------------------------------------------------------
/*
void UpdateEdge(const i3d::Image3d<i3d::GRAY16> &scene,
					 const i3d::Neighbourhood &nbh,
					 std::set<i3d::Vector3d<int> > &edgeSet,
					 i3d::Vector3d<int> &updVoxel,
					 bool periodic)
{
	// for storing values of voxel's neighbourhood
	std::vector<const i3d::GRAY16 *> winNbh;
	// iterate through all neighbours of updVoxel	
	for (size_t i = 1; i < nbh.size(); i++) 
	{
		i3d::Vector3d<int> currVoxel(updVoxel + nbh.offset[i]);

		// handle boundary condition
		if(!ValidateCoords(currVoxel, scene.GetSize(), periodic))
		{
			continue;
		}

		if(scene.GetVoxel(currVoxel) == ID_ECM)
		{
			continue;
		}

		//storing values of currVoxel's neighbourhood
		i3d::GetWindow(scene, 
								currVoxel.x,
								currVoxel.y,
								currVoxel.z,
								nbh, winNbh);

		// (flag = true) -> voxel is "EDGE VOXEL"
		bool isEdgeVoxel = false;

		for (size_t j = 1; j < winNbh.size(); j++) 
		{
			// analyzing if currVoxel is "EDGE VOXEL"
			// compare with neighbourhood values
			// (nbh value NOT ECM and nbh, currVoxel value are DIFFERENT) ->
			//	-> currVoxel is "EDGE VOXEL"
			if ((*winNbh[j] != ID_ECM) && 
				 (scene.GetVoxel(currVoxel) != *winNbh[j])) 
			{
				edgeSet.insert(currVoxel);
				isEdgeVoxel = true;
				//break because we already know that currVoxel is EDGE VOXEL
				break;
			}
		}
		// if voxel doesn't create edge, erase from edgeSet if contains it
		if (!isEdgeVoxel) 
		{

			std::set<i3d::Vector3d<int> >::iterator it=edgeSet.find(currVoxel);

			if(it != edgeSet.end())
			{
				// voxel in edgeSet, we erase it
				edgeSet.erase(it);
			}	
		}

	}
}
*/
//-----------------------------------------------------------------------

