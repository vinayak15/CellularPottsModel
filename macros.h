/**********************************************************************
* 
* macros.h
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
* Description: Collection of useful macros and functions.
*
***********************************************************************/

#ifndef _MACROS_H_
#define _MACROS_H_

#define MAX_STRLEN 256
#define DEG_TO_RAD(x) ((x)*M_PI/180)
#define SQR(x) ((x)*(x))
#define SQRd(x) (((double)(x))*((double)(x)))
#define EVEN(x) (((x) % 2) == 0)
#define ODD(x) (((x) % 2) == 1)

#include <iostream>
#include <string.h>

#define _SHORT_FILE_ strrchr(__FILE__, '/') ? \
			strrchr(__FILE__, '/') + 1 : __FILE__

/// helper macro to unify reports
#define REPORT(x) std::cout \
   	<< std::string(_SHORT_FILE_) << "::" << std::string(__FUNCTION__) \
	<< "(): " << x << std::endl;

/**
 * helper macro to unify debug reports
 *
 * Reports are suppressed in release version of program.
 */

   #define DEBUG_REPORT(x) REPORT(x)

  // #define DEBUG_REPORT(x)


/// helper macro to unify error reports
#define ERROR_REPORT(x) std::string(_SHORT_FILE_)+"::"+__FUNCTION__+"(): "+x

#endif
