/**********************************************************************
* 
* iniparser.cc
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
* Description: A module that can parse extended-INI files. 
* The 'extended' means recognition of couples and triplets as well.
*
***********************************************************************/


#include <iostream>
#include <stdio.h>
#include "iniparser.h"
#include "iniparser.tab.hh"

using namespace std;

extern FILE *yyin;
extern int yyparse (p_section_seq *data);

//---------------------------------------------------------------------------
IniHandler::IniHandler(const char *fname)
{
	 file_name = fname;
#ifdef DEBUG
	 cout << "opening file " << fname << endl;
#endif

	 if ((yyin = fopen(fname, "r")) == NULL)
		  throw "Error: Cannot read file!";

#ifdef DEBUG
	 cout << "file successfully opened" << endl;
	 cout << "parsing file data" << endl;
	 cout << "++++++++++++++++++++++" << endl;
#endif
 
	 p_section_seq data = NULL;

	 if (yyparse(&data) != 0)
		  throw "Error: inifile parse error";

	 *((p_section_seq)this) = *data;

	 delete data;

#ifdef DEBUG
	 cout << "++++++++++++++++++++++" << endl;
	 cout << "data successfully parsed" << endl;
#endif

	 fclose(yyin);

#ifdef DEBUG
	 cout << "file closed" << endl;
	 cout << "++++++++++++++++++++++" << endl;
#endif
}

//---------------------------------------------------------------------------
/*const char* IniHandler::GetFileName() const
{
	 return file_name.c_str();
}*/

//---------------------------------------------------------------------------
IniHandler::~IniHandler() 
{
/*	if (data)
		delete data; */
}

//---------------------------------------------------------------------------
/*t_item_seq& IniHandler::operator[](const char *key)
{
	 t_section_seq::iterator it = data->find(key);

	 // if this section does not exist, throw an exception
	 // segmentation fault
	 if (it == data->end())
	 {
		  string message = string("Error: Cannot find section [") + 
					 key + "]. No such section exists.";
		  throw message.c_str();
	 }

	 // if there is more than one occurrence of 'key' no one
	 // knows which should be taken. Hence throw an exception
	 if (amount(key) > 1)
		  throw "Error: too many identical objects.";

	 return it->second;
}

//---------------------------------------------------------------------------
t_item_seq& IniHandler::operator()(const char *key, size_t order)
{
	 t_section_seq::iterator it = data->lower_bound(key);

	 // if the section with given order does not exist, throw an exception
	 // object to prevent segmentation fault
	 if ((order > (amount(key)-1)) || (it == data->end()))
	 {
		  string message = string("Error: Cannot find section [") + 
					 key + "]. No such section exists.";
		  throw message.c_str();
	 }

	 while (order != 0)
	 {
		  order--;
		  it++;
	 }

	 return (it->second);
}

//---------------------------------------------------------------------------
size_t IniHandler::amount(const char *key)
{
	 return data->count(key);
}

//---------------------------------------------------------------------------
bool IniHandler::find(const char *key)
{
	 return (data->find(key) != data->end());
}
*/

//---------------------------------------------------------------------------

std::ostream& operator<<(std::ostream &os, const t_triplet &t)
{
	os << "(" << t.x << "," << t.y << "," << t.z << ")";
	return os;
}

//---------------------------------------------------------------------------

std::ostream& operator<<(std::ostream &os, const t_couple &t)
{
	os << "(" << t.x << "," << t.y << ")";
	return os;
}

//---------------------------------------------------------------------------


