/**********************************************************************
*
* cpm.cpp
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
* Description: A basic class defining the cellular Potts model in 3D
* extended by diffusion equations and further add-ons to enable
* realistic generation of cellular networks.
*
***********************************************************************/


// Here we consider boundary are always periodic
#ifndef BOUNDARY_PERIODIC
#define BOUNDARY_PERIODIC
#endif



#include <i3d/image3d.h>

#include "cpm.h"
#include "settings.h"
#include "initial.h"
#include "rnd_generators.h"
#include "edge_detection.h"
#include "cell.h"
#include "boundary.h"
#include <valarray>
#include <assert.h>
#include <i3d/regions.h>
#include <i3d/threshold.h>
#include "edge_detection.h"
#include<cuda.h>
#include<cuda_runtime.h>
#include "device_launch_parameters.h"
#include<device_atomic_functions.h>
#include<device_functions.h>
#include<Windows.h>

//-----------------------------------------------------------------------
CPM::CPM(IniHandler *cfg): params(cfg)
{
	 currentStep = 0;
	 sumOfPlannedSteps =
		(int) (*params)["cellular potts model"]["number of steps"];
	 int no_nbh =
		(int) (*params)["cellular potts model"]["neighbourhood"];
	 boltzmannTemp =
		(float) (*params)["cellular potts model"]["boltzmann temperature"];
	 sliceOrder =
		(int) (*params)["rendering"]["order of visualized slice"];

	 switch(no_nbh)
	 {
		  case 6	:
					nbh = i3d::nb3D_o6;
					break;
		  case 18	:
					nbh = i3d::nb3D_o18;
					break;
		  case 26	:
					nbh = i3d::nb3D_o26;
					break;
		  default	:
					DEBUG_REPORT("Neighbourhood set to default "
										"value (6 neighbours) ...");
					nbh = i3d::nb3D_o6;
	 }

	 // Establish ADHESION STRENGTHS
	 CreateTableOfAdhesions(adhesions, params);

	 // Read SHAPE CONSTRAINTS
	 targetVolume =
		(int) (*params)["cellular potts model"]["target volume"];
	 lambdaVolume =
		(int) (*params)["cellular potts model"]["lambda volume"];
	 lambdaSurface =
		(int) (*params)["cellular potts model"]["lambda surface"];

	 // Read relaxation time
	 relaxation = 
		(int) (*params)["cellular potts model"]["relaxation time"];
	 // Read CHEMICAL PARAMETERS
	 chemotaxis = (double) (*params)["pde"]["chemotaxis"];
	 secreteRate = (double) (*params)["pde"]["secr_rate"];
	 decayRate = (double) (*params)["pde"]["decay_rate"];
	 diffCoeff = (double) (*params)["pde"]["diff_coeff"];
	 diffTime = (double) (*params)["pde"]["dt"];
	 diffSpace = (double) (*params)["pde"]["dx"];
	 pdeIters = (int) (*params)["pde"]["pde iterations"];

	 edgeContainer = NULL;
}

//-----------------------------------------------------------------------
CPM::~CPM()
{
	// Deallocate cell structures
	for(unsigned int i = 0; i < cells.size(); i++)
	{
		delete cells[i];
	}

	if (edgeContainer)
		delete edgeContainer;
}

//-----------------------------------------------------------------------
void CPM::ShowCellVolumes() const
{

	 std::cout << "Expected volume of each cell:" << targetVolume << std::endl;
	 for (size_t i=0; i<cells.size(); i++)
	 {
		  std::cout << "(" << i << ":" << cells[i]->Volume() << ")";
	 }
	 std::cout << std::endl;

}

//-----------------------------------------------------------------------
void CPM::InitializePopulation(i3d::Image3d<i3d::GRAY16> &img)
{
	 DEBUG_REPORT("Creating the new population ...");

	 // Initialize an empty scene (3d image)
	 CreateEmptyScene(imgCellIDs, *params);

	 // Fill the scene with the initial number of cells (rough masks)
	 GenerateInitialPopulation(imgCellIDs, *params);

	 // Create cell structure for every cell_id
	 for (size_t i = 0; i <= imgCellIDs.GetMaxValue(); i++)
	 {
		 cells.push_back(new Cell());
	 }
	 MeasureVolumeOfCells();
	 DEBUG_REPORT("Creation of a new population completed.");

	 imgCellIDs.SaveImage("_newinitial.ics");

	 /**
	  * Initialize also the memory buffer for final image rendering.
	  * If we allocated it just before rendering, i.e. many times,
	  * we will slow down the whole simulation process.
	  */
	 imgRendered.CopyMetaData(imgCellIDs);

	 // Initialize chemoatractant's concentration plane
	 imgConcentration.CopyMetaData(imgCellIDs);
	 imgAltConcentration.CopyMetaData(imgCellIDs);
	 img = imgCellIDs;

}

//-----------------------------------------------------------------------
void CPM::ImposeInitialPopulation(i3d::Image3d<i3d::GRAY16> &img)
{
	 imgCellIDs = img;

	 // Create cell structure for every cell_id
	 for (size_t i = 0; i <= imgCellIDs.GetMaxValue(); i++)
	 {
		 cells.push_back(new Cell());
	 }
	 MeasureVolumeOfCells();

	 /**
	  * Initialize also the memory buffer for final image rendering.
	  * If we allocated it just before rendering, i.e. many times,
	  * we will slow down the whole simulation process.
	  */
	 imgRendered.CopyMetaData(imgCellIDs);

	 // Initialize chemoatractant's concentration plane
	 imgConcentration.CopyMetaData(imgCellIDs);
	 imgAltConcentration.CopyMetaData(imgCellIDs);

	 DEBUG_REPORT("Initial cell population read from file you provided.");
}

//-----------------------------------------------------------------------
void CPM::MeasureVolumeOfCells()
{
	// Set volumes of cells
	for (size_t i = 0; i < imgCellIDs.GetImageSize(); i++)									//Error in code
	{
		if (GetObjectType(imgCellIDs.GetVoxel(i)) == TypeCell) 
		{
			cells[imgCellIDs.GetVoxel(i)]->IncrementVolume();
		}
	}
}

//-----------------------------------------------------------------------
void CPM::PrecomputeEdges()
{
	 edgeContainer = new Edges(this->imgCellIDs, this->nbh);
}

//-----------------------------------------------------------------------
void CPM::DoNextStep(float * imgCellIDs, float *imgConcentration)
{
	DEBUG_REPORT("Current step is: " << currentStep);

	size_t rand;

	for (size_t i = 0; i < edgeContainer->size(); i++)
	{
		// 1. Select RANDOM 'source' voxel from edgeSet
		rand = GetRandomUniform(0, edgeContainer->size() - 1);
		//std::set<i3d::Vector3d<int> >::const_iterator it(edgeSet.begin());

		/// OBSOLETE (too slow)
		// 'advance' the iterator RAND times (linear complexity)
		//std::advance(it, rand);

		// random selection of edge voxel with O(1) complexity
		size_t j = edgeContainer->GetIndex(rand);
		i3d::Vector3d<int> sourceCoor = i3d::Vector3d<size_t>(this->imgCellIDs.GetX(j), this->imgCellIDs.GetY(j), this->imgCellIDs.GetZ(j));

		// 2. Load the value of the source voxel
		i3d::GRAY16 sourceValue = imgCellIDs[this->imgCellIDs.GetIndex(sourceCoor)];

		if (sourceValue == ID_ECM)
			DEBUG_REPORT("CHYBA!");

		// 3. Select RANDOM neighbour of selected voxel
		//
		// If we stay in the relaxation time, do not leave the selection
		// process to be random. Prefer the bottom neighbours.
		int selection = GetRandomUniform(1, nbh.size() - 1);

		if ((int)currentStep <= relaxation)
		{
			i3d::Vector3d<int> pos = nbh.offset[selection];

			// The new z-position should be deeper or the same. If not, repeat the 
			// selection process.
			while (pos.z > -1)
			{
				selection = GetRandomUniform(1, nbh.size() - 1);
				pos = nbh.offset[selection];
			}
		}

		i3d::Vector3d<int> targetCoor = sourceCoor + nbh.offset[selection];

		// 4. Check boundary condition
		if (!ValidateCoords(targetCoor, this->imgCellIDs.GetSize()))
		{
			continue;
		}

		// 5. Load the value of the target voxel
		i3d::GRAY16 targetValue = imgCellIDs[this->imgCellIDs.GetIndex(targetCoor)];

		// The value of source voxel is always some cell ID! We need to check
		// the value of the target voxel. It cannot be ID_ECM (it is a solid
		// material and the cells cannot penetrate it). Additionally, it cannot
		// bear the same value as the source voxel does. If so, no change
		// happens.
		if ((targetValue != ID_ECM) && (sourceValue != targetValue))
		{
			if (targetValue == ID_ECM)
				DEBUG_REPORT("Divne targetValue");
			if (sourceValue == ID_ECM)
				DEBUG_REPORT("Divne sourceValue");

			assert((targetValue != ID_ECM) && (sourceValue != ID_ECM));
			double deltaH = 0;

			// propose the spin flip and compute the difference of Hamiltonian
			deltaH = ComputeDeltaH(sourceCoor, targetCoor, imgCellIDs, imgConcentration);

			// is the change accepted?
			if (ProbabilityToSpin(deltaH))
			{
				// do the change
				PerformSpin(sourceCoor, targetCoor, imgCellIDs);

				// update the list of edges
				edgeContainer->Update(this->imgCellIDs, nbh, sourceCoor, imgCellIDs);
			}
		}
	}

	// keep this incrementation at the end of this method!
	currentStep++;
}

//-----------------------------------------------------------------------
bool CPM::ProbabilityToSpin(double deltaH)
{
	double prob;
	//if deltaH = 0 => spin will be performed with 100 percent probability
	if (deltaH <= 0)
	{
		return true;
	}
	else
	{
		//computing probability of spin based on deltaH
		//larger deltaH means larger probability
		prob = exp(-(deltaH / boltzmannTemp));
		return GetRandomUniform(0.0f, 1.0f) < prob;
	}
}

//-----------------------------------------------------------------------
void CPM::PerformSpin(i3d::Vector3d<int> source, i3d::Vector3d<int> target, float * imgCellIDs)
{
	i3d::GRAY16 sourceID = imgCellIDs[this->imgCellIDs.GetIndex(source)];
	i3d::GRAY16 targetID = imgCellIDs[this->imgCellIDs.GetIndex(target)];

	cells[sourceID]->DecrementVolume();
	cells[targetID]->IncrementVolume();

	imgCellIDs[this->imgCellIDs.GetIndex(source)] = imgCellIDs[this->imgCellIDs.GetIndex(target)];
}

//-----------------------------------------------------------------------
int CPM::LocalSurfaceAfterChange(i3d::Vector3d<int> changingCoor,
	int cellID, int newID, float *imgCellIDs)
{
	//	localSurface -> number of voxels with cellID that create surface 
	//					after performing change
	//				 -> these voxels belong to neighbourhood of changingCoor 
	//					(changingCoor included)
	int localSurface = 0;
	std::vector<const float *> winNbh;
	//	iterate through all neighbour voxels of changingCoor voxel
	for (size_t i = 0; i < nbh.size(); i++) {
		i3d::Vector3d<int> nbhCoor(changingCoor + nbh.offset[i]);
		// handle boundary conditions
		if (!ValidateCoords(nbhCoor, this->imgCellIDs.GetSize()))
		{
			continue;
		}
		//	skip current voxels with another ID than cellID
		if (imgCellIDs[this->imgCellIDs.GetIndex(nbhCoor)] != cellID)
		{
			continue;
		}
		// This condition is here because of OPTIMIZATION.
		// At this point we know that "current voxel" has id=cellID and
		// is neighbour of voxel with coordinates=changeCoor.
		// ChangeCoor voxel has id=newID. If newID is different from cellID, 
		// we can say "current voxel"(the one with id=cellID) belongs to SURFACE
		if (imgCellIDs[this->imgCellIDs.GetIndex(nbhCoor)] != newID)
		{
			localSurface++;
			continue;
		}
		// If program gets to this place in code, 
		// it means newID = cellID => current voxel - cellID
		//							  changeCoor voxel - cellID
		// We then have to check whole neighbourhood of current voxel and
		// decide if it creates SURFACE
		//take neighbourhood of site

		i3d::Neighbourhood rnb;
		size_t tmp = GetNbh(this->imgCellIDs, nbhCoor.x, nbhCoor.y, nbhCoor.z, nbh, rnb);
		winNbh.resize(nbh.size());
		i3d::VectContainer::const_iterator off;
		int q = 0;
		for (off = nbh.offset.begin(); off != nbh.offset.end(); ++off)
		{
			winNbh[q++] = &imgCellIDs[nbhCoor.x + off->x + (nbhCoor.y + off->y)*this->imgCellIDs.GetSizeX() + (nbhCoor.z + off->z)*this->imgCellIDs.GetSizeX() * this->imgCellIDs.GetSizeY()];
		}


		for (size_t j = 1; j < winNbh.size(); j++)
		{
			// check if site is on surface
			if (imgCellIDs[this->imgCellIDs.GetIndex(nbhCoor)] != *winNbh[j])
			{
				localSurface++;
				break;
			}
		}
	}
	return localSurface;
}

//-----------------------------------------------------------------------
double CPM::ComputeDeltaH(i3d::Vector3d<int> source, i3d::Vector3d<int> target, float *imgCellIDs, float * imgConcentration)
{
	double H_before = 0.0;
	double H_after = 0.0;

	i3d::GRAY16 sourceID = imgCellIDs[this->imgCellIDs.GetIndex(source)];
	i3d::GRAY16 targetID = imgCellIDs[this->imgCellIDs.GetIndex(target)];
	i3d::GRAY8 nbhID;

	// *** H_adhesion ***
	//
	// In the following loop we need to inspect, how the possible change
	// affects the relations ship of all the neighbouring voxels.
	//
	for (size_t i = 1; i < nbh.size(); i++)
	{
		i3d::Vector3d<int> nbhCoor = source + nbh.offset[i];
		// handle boundary condition
		if (!ValidateCoords(nbhCoor, this->imgCellIDs.GetSize()))
		{
			continue;
		}
		nbhID = imgCellIDs[this->imgCellIDs.GetIndex(nbhCoor)];

		H_after += adhesions[GetObjectType(targetID)][GetObjectType(nbhID)];
		H_before += adhesions[GetObjectType(sourceID)][GetObjectType(nbhID)];
	}

	// *** H_shape ***
	//
	// VOLUME CONSTRAINT
	//
	// Here, we inspect, how the spin spin affects the volume of the object
	// which the source voxel belongs to.
	//
	if (targetID == ID_MEDIUM)
	{
		//situation -> source CELL, target MEDIUM
		//if change happens, volume of cell will decrease by one
		H_before += lambdaVolume *
			SQRd((cells[sourceID]->Volume() - TargetVolume()));

		H_after += lambdaVolume *
			SQRd((cells[sourceID]->Volume() - 1 - TargetVolume()));

	}
	else if (sourceID == ID_MEDIUM)
	{
		//situation -> source MEDIUM, target CELL
		//if change happens, volume of cell will increase by one
		H_before += lambdaVolume *
			SQRd((cells[targetID]->Volume() - TargetVolume()));

		H_after += lambdaVolume *
			SQRd((cells[targetID]->Volume() + 1 - TargetVolume()));
	}
	else
	{
		//situation -> source CELL1, target CELL2
		//if change happens, volume of cell1 will decrease by one
		//					 volume of cell2 will increase by one
		H_before += lambdaVolume * (
			SQRd((cells[sourceID]->Volume() - TargetVolume())) +
			SQRd((cells[targetID]->Volume() - TargetVolume())));
		H_after += lambdaVolume * (
			SQRd((cells[sourceID]->Volume() - 1 - TargetVolume())) +
			SQRd((cells[targetID]->Volume() + 1 - TargetVolume())));

	}

	//
	// SURFACE CONSTRAINT
	// 
	if (targetID == ID_MEDIUM)
	{
		// get local surface of cell with id=sourceID before performing spin 
		// change on coordinate "source"
		int deltaSurBefore = LocalSurfaceAfterChange(source, sourceID, sourceID, imgCellIDs);
		// get local surface of cell with id=sourceID after performing spin 
		// change on coordinate "source"
		int deltaSurAfter = LocalSurfaceAfterChange(source, sourceID, targetID, imgCellIDs);
		//situation -> source CELL, target MEDIUM
		H_before += lambdaSurface * deltaSurBefore;
		H_after += lambdaSurface * deltaSurAfter;
	}
	else if (sourceID == ID_MEDIUM)
	{
		// get local surface of cell with id=targetID before performing spin 
		// change on coordinate "source"	
		int deltaSurBefore = LocalSurfaceAfterChange(source, targetID, sourceID, imgCellIDs);
		// get local surface of cell with id=targetID after performing spin 
		// change on coordinate "source"
		int deltaSurAfter = LocalSurfaceAfterChange(source, targetID, targetID, imgCellIDs);
		//situation -> source CELL, target MEDIUM
		H_before += lambdaSurface * deltaSurBefore;
		H_after += lambdaSurface * deltaSurAfter;
	}
	else
	{
		// get local surface of cell with id=sourceID on coordinate "source" and
		// local surface of cell with id=targetID on coordinate "target" BEFORE
		// performing spin change
		int deltaSurBefore = LocalSurfaceAfterChange(source, sourceID,
			sourceID, imgCellIDs) + LocalSurfaceAfterChange(target, targetID, targetID, imgCellIDs);
		// get local surface of cell with id=sourceID on coordinate "source" and
		// local surface of cell with id=targetID on coordinate "target" AFTER
		// performing spin change
		int deltaSurAfter = LocalSurfaceAfterChange(source, sourceID, targetID, imgCellIDs);
		// to get local surface of cell with id=targetID after spin change, 
		// we need to change voxel value on coordinate "source" to targetID
		// and then after applying function LocalSurfaceAfterChange return it
		// back again
		imgCellIDs[this->imgCellIDs.GetIndex(source)] = targetID;
		deltaSurAfter += LocalSurfaceAfterChange(target, targetID, targetID, imgCellIDs);
		imgCellIDs[this->imgCellIDs.GetIndex(source)] = sourceID;
		//situation -> source CELL1, target CELL2
		H_before += lambdaSurface * deltaSurBefore;
		H_after += lambdaSurface * deltaSurAfter;
	}

	/// *** H_chemical ***
	double DDH = 0;

	DDH = (chemotaxis * (imgConcentration[this->imgCellIDs.GetIndex(source)] -
		imgConcentration[this->imgCellIDs.GetIndex(target)]));

	return (H_after - H_before) - DDH;
}

//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void CPM::Secrete() 
{
	const double increase = secreteRate * diffTime;
	const double decay = (1.0f - decayRate);
	const size_t sz = imgCellIDs.GetImageSize();
	const i3d::GRAY16 *currIDPtr = imgCellIDs.GetFirstVoxelAddr();
	float *currConcPtr = imgConcentration.GetFirstVoxelAddr();

	for (size_t i=0; i<sz; i++)
	{
		 // Only the cells secrete chemoattractant, medium and ECM does not
		 if ((*(currIDPtr) != ID_MEDIUM) && (*(currIDPtr) != ID_ECM))
		 {
			  *(currConcPtr) += increase;
		 }
		 else // Outside the cells, chemoattractant decays
		 {
			  *(currConcPtr) *= decay;
		 }

		 currIDPtr++;
		 currConcPtr++;
	}

}

//---------------------------------------------------------------------------------------------


__global__ void CudaSecrete(float *imgCellIDs, float *imgConcentration, float secreteRate, float diffTime, float decayRate)      //Decalred global so that it can access by any class
{
	size_t index = threadIdx.x + blockIdx.x*blockDim.x;																//threadIdx.x is index of thread in x direction and blockIdx.x is index of block in X direction  and blockDim.x is dimension of block in x direction
	const double increase = secreteRate*diffTime;							
	const double decay=(1.0f - decayRate);
	if (imgCellIDs[index] != ID_MEDIUM && imgCellIDs[index] != ID_ECM)			
		imgConcentration[index] += increase;
	else
		imgConcentration[index] *= decay;
}

//---------------------------------------------------------------------------------------------
__global__ void CudaDiffuse(float *imgConcentration, float *imgAltConcentration, int width, int height, int depth, int sliceSize, double diffConstant)
{
	size_t index = threadIdx.x + blockIdx.x*blockDim.x;
	size_t x = index%width;
	size_t y = (index / width) % height;
	size_t z = ((index / width) / height);
	size_t left, right, above, below, front, behind;
	float sum;
	left = index - 1;																//Setting left,right, above, below,front, behind index for each and every thread.
	right = index + 1;
	above = index - width;																
	below = index + width;
	front = index - sliceSize;
	behind = index + sliceSize;

#ifdef BOUNDARY_PERIODIC															//Treating each and every boundary condition seprately
	if ( x == 0 &&  y == (height - 1) &&  z<(depth - 1) &&  z>0) {		//here block Idx.x is index fof image in x direction  
		left = index + (width - 1);
		below = index - width*(height - 1);															// y is index of image in y direction
		// diffuse chemoattractant in every voxel (except for boundary ones)						// z is index of image in z direction
		sum = (-6.0f) * (imgConcentration[index]);
		sum += (imgConcentration[left]) + (imgConcentration[right]) + (imgConcentration[above]) + (imgConcentration[below]) + (imgConcentration[front]) + (imgConcentration[behind]);

		// store the result
		imgAltConcentration[index] = imgConcentration[index] + sum * diffConstant;
	}
	else if ( x == (width - 1) &&  (y == height - 1) &&  z<depth - 1 &&  z>0) {
		right = index - (width - 1);
		below = index - width*(height - 1);

		// diffuse chemoattractant in every voxel (except for boundary ones)
		sum = (-6.0f) * (imgConcentration[index]);
		sum += (imgConcentration[left]) + (imgConcentration[right]) + (imgConcentration[above]) + (imgConcentration[below]) + (imgConcentration[front]) + (imgConcentration[behind]);

		// store the result
		imgAltConcentration[index] = imgConcentration[index] + sum * diffConstant;
	}
	else if ( x == (width - 1) &&  y == 0 &&  z<depth - 1 &&  z>0) {
		right = index - (width - 1);
		above = index + width*(height - 1);
		// diffuse chemoattractant in every voxel (except for boundary ones)
		sum = (-6.0f) * (imgConcentration[index]);
		sum += (imgConcentration[left]) + (imgConcentration[right]) + (imgConcentration[above]) + (imgConcentration[below]) + (imgConcentration[front]) + (imgConcentration[behind]);

		// store the result
		imgAltConcentration[index] = imgConcentration[index] + sum * diffConstant;
	}
	else if ( x == 0 &&  y == 0 &&  z<(depth - 1) &&  z>0) {
		left = index + (width - 1);
		above = index + width*(height - 1);
		// diffuse chemoattractant in every voxel (except for boundary ones)
		sum = (-6.0f) * (imgConcentration[index]);
		sum += (imgConcentration[left]) + (imgConcentration[right]) + (imgConcentration[above]) + (imgConcentration[below]) + (imgConcentration[front]) + (imgConcentration[behind]);

		// store the result
		imgAltConcentration[index] = imgConcentration[index] + sum * diffConstant;
	}
	else if ( x == (width - 1) &&  z<(depth - 1) &&  z>0) {
		right = index - (width - 1);
		// diffuse chemoattractant in every voxel (except for boundary ones)
		sum = (-6.0f) * (imgConcentration[index]);
		sum += (imgConcentration[left]) + (imgConcentration[right]) + (imgConcentration[above]) + (imgConcentration[below]) + (imgConcentration[front]) + (imgConcentration[behind]);

		// store the result
		imgAltConcentration[index] = imgConcentration[index] + sum * diffConstant;
	}
	else if ( x == 0 &&  z<depth - 1 &&  z>0) {
		left = index + (width - 1);
		// diffuse chemoattractant in every voxel (except for boundary ones)
		sum = (-6.0f) * (imgConcentration[index]);
		sum += (imgConcentration[left]) + (imgConcentration[right]) + (imgConcentration[above]) + (imgConcentration[below]) + (imgConcentration[front]) + (imgConcentration[behind]);

		// store the result
		imgAltConcentration[index] = imgConcentration[index] + sum * diffConstant;
	}
	else if ( y == height - 1 &&  z<depth - 1 &&  z>0) {
		below = index - width*(height - 1);
		// diffuse chemoattractant in every voxel (except for boundary ones)
		sum = (-6.0f) * (imgConcentration[index]);
		sum += (imgConcentration[left]) + (imgConcentration[right]) + (imgConcentration[above]) + (imgConcentration[below]) + (imgConcentration[front]) + (imgConcentration[behind]);

		// store the result
		imgAltConcentration[index] = imgConcentration[index] + sum * diffConstant;
	}
	else if ( y == 0 &&  z<depth - 1 &&  z>0) {
		above = index + width*(height - 1);
		// diffuse chemoattractant in every voxel (except for boundary ones)
		sum = (-6.0f) * (imgConcentration[index]);
		sum += (imgConcentration[left]) + (imgConcentration[right]) + (imgConcentration[above]) + (imgConcentration[below]) + (imgConcentration[front]) + (imgConcentration[behind]);

		// store the result
		imgAltConcentration[index] = imgConcentration[index] + sum * diffConstant;
	}
	else if ( (x >= 1) && ( y >= 1) && ( x<(width - 1)) &&  (y<(height - 1)) &&  (z<(depth - 1)) &&  (z>0) )
	{
		// diffuse chemoattractant in every voxel (except for boundary ones)
		sum = (-6.0f) * (imgConcentration[index]);
		sum += (imgConcentration[left]) + (imgConcentration[right]) + (imgConcentration[above]) + (imgConcentration[below]) + (imgConcentration[front]) + (imgConcentration[behind]);

		// store the result
		imgAltConcentration[index] = imgConcentration[index] + sum * diffConstant;
	}
	#endif
												
	imgConcentration[index] = imgAltConcentration[index];
	
}


//-----------------------------------------------------------------------
void CPM::Render(float *imgCellIDs)
{
	 DEBUG_REPORT("Creating the output image");
	 for (size_t i = 0; i < this->imgCellIDs.GetImageSize(); i++)
		 this->imgCellIDs.SetVoxel(i, imgCellIDs[i]);
	 const char *colors = (*params)["rendering"]["true colors"];

	 if (strcmp(colors, "true") == 0)
	 {
		GrayToRGB(this->imgCellIDs, this->imgCellIDs, this->imgCellIDs, imgRendered);
	 }
	 else if (strcmp(colors, "false") == 0)
	 {
		  i3d::RGB16 rgbValue;

		  for (size_t i=0; i<imgRendered.GetImageSize(); i++)
		  {
				i3d::GRAY16 cellID = this->imgCellIDs.GetVoxel(i);

				if (cellID == ID_MEDIUM)
				{
					 rgbValue.red = rgbValue.green = rgbValue.blue = 0;
				}
				else if (cellID == ID_ECM)
				{
					 rgbValue.red = rgbValue.green = rgbValue.blue = 128;
				}
				else // the places where the cells are located
				{
					 // yellow color
					 rgbValue.red = 255;
					 rgbValue.green = 255;
					 rgbValue.blue = 0;
				}

				imgRendered.SetVoxel(i, rgbValue);
		  }
	 }
	 else
	 {
		  throw ERROR_REPORT("Unknown key value. Expected 'true'/'false'.");
	 }


	 DEBUG_REPORT("Image completed");
}

void CPM::Render()
{
	DEBUG_REPORT("Creating the output image");

	const char *colors = (*params)["rendering"]["true colors"];

	if (strcmp(colors, "true") == 0)
	{
		GrayToRGB(imgCellIDs, imgCellIDs, imgCellIDs, imgRendered);
	}
	else if (strcmp(colors, "false") == 0)
	{
		i3d::RGB16 rgbValue;

		for (size_t i = 0; i<imgRendered.GetImageSize(); i++)
		{
			i3d::GRAY16 cellID = imgCellIDs.GetVoxel(i);

			if (cellID == ID_MEDIUM)
			{
				rgbValue.red = rgbValue.green = rgbValue.blue = 0;
			}
			else if (cellID == ID_ECM)
			{
				rgbValue.red = rgbValue.green = rgbValue.blue = 128;
			}
			else // the places where the cells are located
			{
				// yellow color
				rgbValue.red = 255;
				rgbValue.green = 255;
				rgbValue.blue = 0;
			}

			imgRendered.SetVoxel(i, rgbValue);
		}
	}
	else
	{
		throw ERROR_REPORT("Unknown key value. Expected 'true'/'false'.");
	}


	DEBUG_REPORT("Image completed");
}


//-----------------------------------------------------------------------
void CPM::StoreToFile() const
{
	 const char *agreement = (*params)["rendering"]["store data to disk"];

	 if (strcmp(agreement,"true") == 0)
	 {
		  char fname[MAX_STRLEN];
		  sprintf(fname, "img_%.4lu.ics", currentStep);

		  DEBUG_REPORT("Saving the image file: " << fname);
		  imgRendered.SaveImage(fname);
		  DEBUG_REPORT("File saved successfully.");

		  /*
		  // TODO: smazat po odladeni
		  // begin - koncentrace (separe) 
		  char fname2[MAX_STRLEN];
		  sprintf(fname2, "signals_%.4lu.ics", currentStep);
		  DEBUG_REPORT("Saving the image file: " << fname2);
		  imgConcentration.SaveImage(fname2);
		  DEBUG_REPORT("File saved successfully.");
		  // end - koncentrace
		  */
	 }
	 else if (strcmp(agreement,"false") == 0)
	 {
		  // do nothing
	 }
	 else
	 {
		  throw ERROR_REPORT("Unknown key value. Expected 'true'/'false'");
	 }
}
//-----------------------------------------------------------------------
void CPM::GetImage(float *img) {
	for (size_t i = 0; i < imgCellIDs.GetImageSize(); i++)
		img[i] = imgCellIDs.GetVoxel(i);
}

//-----------------------------------------------------------------------

void CPM::SetImage(float *imgCon) {
	for (size_t i = 0; i < imgCellIDs.GetImageSize(); i++) {
		imgConcentration.SetVoxel(i, imgCon[i]);
	}

}

//-----------------------------------------------------------------------
ObjectType GetObjectType(int objectID)
{
	 if (objectID < 0)
	 {
		  ERROR_REPORT("Invalid ID value!");
	 }

	 switch (objectID)
	 {
		  case ID_MEDIUM: return TypeMedium;
		  case ID_ECM: return TypeECM;
		  default: {};
	 }

	 return TypeCell;
}


//-----------------------------------------------------------------------
