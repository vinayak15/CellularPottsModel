/**********************************************************************
* 
* initial.cpp
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
* Description: Generation of initial cell population. Based on the paper:
*
* "Svoboda D., Ulman V. Towards a Realistic Distribution of Cells 
* in Synthetically Generated 3D Cell Populations. In Image Analysis 
* and Processing - ICIAP 2013, Part II, LNCS 8157.  Berlin, Heidelberg: 
* Springer-Verlag,  pp 429-438, September 2013, ISBN 978-3-642-41183-0"
*
***********************************************************************/

#include <i3d/threshold.h>
#include <i3d/regions.h>
#include <i3d/neighbours.h>

#include "shape_rendering.h"
#include "collisions.h"

#include "initial.h"

#include <fstream>
#include <iomanip>


//---------------------------------------------------------------------------
template <class T>
void CreateEmptyScene(i3d::Image3d<T> &scene, IniHandler &ini)
{
	 t_triplet fsize = ini["specimen"]["size in microns"];
	 t_triplet res = ini["resolution"]["pixels per micron"];

	 i3d::Vector3d<size_t> sz;

	 sz.x = (size_t) floor(fsize.x * res.x);
	 sz.y = (size_t) floor(fsize.y * res.y);
	 sz.z = (size_t) floor(fsize.z * res.z);

	 scene.MakeRoom(sz);
	 scene.SetResolution(i3d::Vector3d<float>(res.x, res.y, res.z));
	 scene.SetOffset(0.0f);

	 for (size_t i=0; i<scene.GetImageSize(); i++)
	 {
		  scene.SetVoxel(i, ID_MEDIUM);
	 }

	 DEBUG_REPORT("The size of scene is going to be: " << 
					  scene.GetSize() << " pixels");
}

//---------------------------------------------------------------------------
template <class T>
void GenerateInitialPopulation(i3d::Image3d<T> &scene, IniHandler &ini)
{
	 int numOfCells = ini["cell population"]["number of cells"];
	 float cellDiameter = ini["cell population"]["mean cell diameter in microns"];

	 i3d::Image3d<bool> occupiedSpace;
	 occupiedSpace.CopyMetaData(scene);

	 DEBUG_REPORT("Creating ECM ...");
	 float thickness = ini["specimen"]["ECM thickness in microns"];
	 size_t height = (size_t) (thickness * scene.GetResolution().GetZ());

	 size_t skip = height*scene.GetSliceSize();
	 for (size_t i=scene.GetImageSize()-skip; i<scene.GetImageSize(); i++)
	 {
		  scene.SetVoxel(i, ID_ECM);
	 }
	 DEBUG_REPORT("ECM completed.");


	 DEBUG_REPORT("Creating individual cell masks ...");

	 for (int i=0; i<numOfCells; i++)
	 {
		  DEBUG_REPORT("Cell order: " << i);

		  i3d::Image3d<bool> cellMask;
		  RenderBlankCellMask(scene.GetResolution().GetRes(), 
									 cellMask, cellDiameter); 

		  i3d::Threshold(scene, occupiedSpace, 
							  T(1), std::numeric_limits<T>::max());
		  
		  i3d::Vector3d<float> *fPos = FindPosition(cellMask, occupiedSpace); 
		  
		  if (fPos == NULL)
		  {
				throw ERROR_REPORT("VOI too small for new cell"); 
		  }

		  i3d::Vector3d<float> 
					 centre = i3d::PixelsToMicrons(cellMask.GetSize(),
															 cellMask.GetResolution())/2.0f;
		  
		  cellMask.SetOffset(*fPos - centre);
		  delete fPos;

		  AddNewLabel(scene, cellMask);
	 }

	 DEBUG_REPORT("Creation of cell masks completed.");


	 scene.SaveImage("_initial_.ics");

}

//---------------------------------------------------------------------------

void CreateTableOfAdhesions(int table[][3], IniHandler *cfg)
{
	DEBUG_REPORT("Creating table of adhesions");

	table[TypeCell][TypeMedium] = 
			  table[TypeMedium][TypeCell] =
			  (*cfg)["adhesion"]["cell-medium"];

	table[TypeCell][TypeECM] = 
			  table[TypeECM][TypeCell] =
			  (*cfg)["adhesion"]["cell-ecm"];

	table[TypeCell][TypeCell] = (*cfg)["adhesion"]["cell-cell"];

	// Zero values by definition as ECM and medium do not move.
	// These values are never used.
	table[TypeMedium][TypeECM] = table[TypeECM][TypeMedium] = 0;
	table[TypeECM][TypeECM] = 0;
	table[TypeMedium][TypeMedium] = 0;


	std::cout << "Printing table of adhesions: " << std::endl;
	for (size_t i=0; i<3; i++)
	{
		 for (size_t j=0; j<3; j++)
		 {
			  std::cout << std::setw(4) << table[j][i];
		 }
		 std::cout << "\n";
	}


}

//---------------------------------------------------------------------------
// explicit instantiation
//---------------------------------------------------------------------------

template void CreateEmptyScene(i3d::Image3d<i3d::GRAY16> &scene, IniHandler &ini);

template void GenerateInitialPopulation(i3d::Image3d<i3d::GRAY16> &scene, IniHandler &ini);
