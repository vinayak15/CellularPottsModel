/**********************************************************************
* 
* initial.h
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

#ifndef _INITIAL_H_
#define _INITIAL_H_

#include <i3d/image3d.h>
#include "ini/iniparser.h"

/**
 * This function loads the data from ini-structure and creates the basic
 * empty scene (labeled image). The result is stored in the first parameter
 * given as a reference.
 */
template <class T>
void CreateEmptyScene(i3d::Image3d<T> &scene, IniHandler &ini);

/**
 * Given the initial empty scene and an ini-structure, this function fills
 * the scene with the expected number of cell. This way, the initial cell
 * population is established.
 */
template <class T>
void GenerateInitialPopulation(i3d::Image3d<T> &scene, IniHandler &ini);

/**
 * Given the initial cell distribution in the image 'cells', the individual
 * cell components are generated into the second image and stored there.
 */
template <class S, class T>
void GenerateCellCompartments(i3d::Image3d<S> &components,
										const i3d::Image3d<T> &cells,
										const i3d::Neighbourhood &nbh);

/**
 * Create a table of static binding energies
 * 1st line -> number of types (all components + medium, ECM)
 * Next lines -> diagonal matrix, starting with 1 element (0 0)ending with
 * n elements
 */
void CreateTableOfAdhesions(int table[][3], IniHandler *cfg);

#endif


