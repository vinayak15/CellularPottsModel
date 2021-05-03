/**********************************************************************
* 
* settings.h
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
* Description: A toolkit configuration file. 
*
***********************************************************************/

#ifndef _SETTINGS_H_
#define _SETTINGS_H_

#include "macros.h"

// basic object types
typedef enum {
	 TypeMedium = 0,
	 TypeECM = 1,
	 TypeCell = 2
} ObjectType;


// object IDs (cells obtain IDs: 2,3,4,...)
#define ID_MEDIUM		0
#define ID_ECM			1

#endif
