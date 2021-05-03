/**********************************************************************
*
* cpm.h
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
* Description: A basic class defining the cellular Potts model in 3D
* extended by diffusion equations and further add-ons to enable
* realistic generation of cellular networks.
*
***********************************************************************/

#ifndef _CPM_H_
#define _CPM_H_

#include <i3d/image3d.h>
#include <i3d/neighbours.h>
#include "ini/iniparser.h"

#include <vector>
#include <set>

#include "settings.h"
#include "cell.h"
#include "edge_detection.h"
#include "edge_detection.h"
#include<cuda.h>
#include<cuda_runtime.h>
#include "device_launch_parameters.h"
#include<device_atomic_functions.h>
#include<device_functions.h>
#include<Windows.h>
/**
 * A function associating the type of object with voxel spin/ID.
 */
ObjectType GetObjectType(int objectID);
__global__ void CudaSecrete(float *, float *, float, float, float);
__global__ void CudaDiffuse(float *, float *, int, int, int,int, double);
/**
 * A main class defining the cellular Potts model
 */
class CPM
{
  public:
			 /// constructor
			 CPM(IniHandler *cfg);

			 /// destructor
			 ~CPM();

			 /**
			  * create an initial cell population based on the
			  * setting given in the configuration file
			  */
			 void InitializePopulation(i3d::Image3d<i3d::GRAY16> &img);

			 /**
			  * Do not generate the new cell population (it is time consuming)
			  * but accept the imposed labeled image file
			  */
			 void ImposeInitialPopulation(i3d::Image3d<i3d::GRAY16> &img);
			 
			 // Set initial volume for all cells
			 void MeasureVolumeOfCells();
			 
			 // Set initial surface for all cells
			 void MeasureSurfaceOfCells();
			 
			 // Get number of cell voxels with spin "cellID" that create surface
			 //	after changing spin on coordinate "changingCoor" to spin "newID"
			 //	Function returns only local number of surface voxels, voxels
			 //	have to be in neighbourhood of changeCoor voxel		 
			 int LocalSurfaceAfterChange(i3d::Vector3d<int> changingCoor,
										int cellID, int newID , float *);

			 /// get the expected length of the simulation
			 size_t GetOverallDuration() const { return sumOfPlannedSteps; };

			 /// get the order of the current simulation step
			 size_t GetCurrentStep() const { return currentStep; };

			 /// perform next simulation step and increase the counter
			 void DoNextStep(float * , float *);

			 /// show the cel population
			 const i3d::Image3d<i3d::GRAY16>& ShowPopulation() const
			 { return imgCellIDs; };

			 /// show the used neighbourhood
			 const i3d::Neighbourhood& GetNeighbourhood() const
			 { return nbh; };

			 /// show the final (appropriately labeled) image
			 const i3d::Image3d<i3d::RGB16>& GetRenderedImageData() const
			 { return imgRendered; };

			 /// render the current configuration into auxiliary memory buffer
			 void Render(float *);
			 void Render();
			 /// save the rendered data into the file
			 void StoreToFile() const;

			 /// get order of slice that should be visualized
			 int GetSliceOrder() const
			 { return sliceOrder; };
			 
			 /// Return cell's target volume
			 int TargetVolume() const {
			 	return targetVolume;
			 }

			 /// Show list of cells with their volumes
			 void ShowCellVolumes() const;
			 
			 ///PDE
			 /// Diffuse chemoatractant to environment
			 void Diffuse();

			 /// Cells secrete chemoatractant to environment
			 void Secrete();
			
			 /// Return number of pde iterations
			 int PdeIterations() const {
			 	return pdeIters;
			 }

			 /// Get relaxation time
			 int GetRelaxationTime() const { return relaxation; }

			 /// Initially compute the edges
			 void PrecomputeEdges();

			 /// Show model parameters
			 const IniHandler *ShowParams() const { return params; }

			 void SetImage(float *);
			 void GetImage(float*);

  private:
			 /// compute difference between Hamiltonian AFTER and BEFORE spin
			 double ComputeDeltaH(i3d::Vector3d<int> source,
								i3d::Vector3d<int> target , float * , float *);

			 /// based on provided deltaH decides if spin will be performed
			 bool ProbabilityToSpin(double deltaH);

			 /// perform spin (assign source's cell ID to target)
			 void PerformSpin(i3d::Vector3d<int> source,
								i3d::Vector3d<int> target , float *);

			 /// return binding energy between two voxels
			 // int BindingEnergy(i3d::GRAY16 cellID1, i3d::GRAY16 cellID2);
			 	 

			 /// order of the current step (the first one is '0')
			 size_t currentStep;

			 /// order of slice that should be visualized in QT frame
			 int sliceOrder;

			 /// relaxation time
			 int relaxation;

			 /// the total expected duration of the simulation
			 size_t sumOfPlannedSteps;
			 
			 ///Boltzmann temperature for CPM system
			 float boltzmannTemp;

			 /// 3D data storing cell ids
			 i3d::Image3d<i3d::GRAY16> imgCellIDs;
			 
			 /// cell's attributes of all cells
			 std::vector<Cell*> cells;

			 /// input parameter given by the user
			 IniHandler *params;
			 
			 /// the rendered image data
			 i3d::Image3d<i3d::RGB16> imgRendered;
 
			 /// detected edge voxels in the scene
			 Edges *edgeContainer;

			 /// neighbourhood used in the system
			 i3d::Neighbourhood nbh;

			 /// adhesion strength between different types of voxels (objects)
			 int adhesions[3][3];
			 
			 /// SHAPE CONSTRAINTS
			 /// target volume of cell
			 int targetVolume;
			 /// factor of volume importance
			 int lambdaVolume;
			 /// factor of volume importance
			 int lambdaSurface;
			 
			 ///CHEMICAL parameters
			 // plane containing concentration of chemoatractant on every voxel
			 i3d::Image3d<float> imgConcentration;
			 
			 // Used as temporary memory in the diffusion step
			 // (addresses will be swapped for every time step, so
			 // never directly use them!!! Access is guaranteed to be correct
			 // through user interface)
			 i3d::Image3d<float> imgAltConcentration;
			
			 double chemotaxis;
			 double secreteRate;
			 double decayRate;

			 // diffusion parameters
			 double diffCoeff;
			 double diffTime;
			 double diffSpace;
			 
			 // Number of pde iterations
			 // (pde consists of SECRETION and DIFFUSION of chemoatractant)
			 double pdeIters;
};



#endif
