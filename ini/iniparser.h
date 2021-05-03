/**********************************************************************
* 
* iniparser.h
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

#ifndef _INIPARSER_
#define _INIPARSER_

#include <stdlib.h>
#include <string.h>
#include <string>
#include <map>
#include <iostream>

#define MAX_STR_LEN 256

//---------------------------------------------------------------------------

typedef struct
	{
		 float x, y, z;
	} t_triplet;

typedef struct 
	{
		 float x, y;
	} t_couple;

typedef union
{
	 char str[MAX_STR_LEN];
	 float f_number;
	 int	i_number;
	 t_triplet triplet;
	 t_couple couple;
} selection;


//---------------------------------------------------------------------------
std::ostream& operator<<(std::ostream &os, const t_triplet &t);
std::ostream& operator<<(std::ostream &os, const t_couple &t);
//---------------------------------------------------------------------------

typedef struct _t_store
{
	 // data store
	 selection value;

	 // type checker - prevents from null pointer assignment
	 enum { my_string, my_real, my_integer, my_couple, my_triplet } my_type;

	 // input methods
	 _t_store& operator=(int i) 
	 { 
		  my_type = my_integer; 
		  value.i_number = i; 
		  return *this; 
	 }

	 _t_store& operator=(float f) 
	 { 
		  my_type = my_real; 
		  value.f_number = f; 
		  return *this; 
	 }

	 _t_store& operator=(t_triplet &t) 
	 { 
		  my_type = my_triplet; 
		  value.triplet = t; 
		  return *this; 
	 }

	 _t_store& operator=(t_couple &t) 
	 { 
		  my_type = my_couple; 
		  value.couple = t; 
		  return *this; 
	 }

	 _t_store& operator=(const char *s) 
	 { 
		  my_type = my_string; 
		  strcpy(value.str, s); 
		  return *this; 
	 }

	 // output methods
	 operator int() const 
	 {
		  if (my_type != my_integer)
		  {
				if (my_type == my_real)
					 return (int)value.f_number;
				else
					 throw "Invalid conversion";
		  }
			
		  return value.i_number; 
	 }

	 operator float() const 
	 {
		  if (my_type != my_real)
		  {
				if (my_type == my_integer)
					 return (float)value.i_number;
				else
					 throw "Invalid conversion";
		  }
		  
		  return value.f_number; 
	 }

	 operator t_triplet() const 
	 {
		  if (my_type != my_triplet)
				throw "Invalid conversion";
		  
		  return value.triplet; 
	 }

	 operator t_couple() const 
	 {
		  if (my_type != my_couple)
				throw "Invalid conversion";
		  return value.couple; 
	 }

	 operator const char*() const 
	 { 
		  if (my_type != my_string)
				throw "Invalid conversion";
		  
		  return value.str; 
	 }

} t_store;

//---------------------------------------------------------------------------

template <class TKey, class TData, class TCompare>
class my_multimap : public std::multimap<TKey, TData, TCompare>
{ 
  typedef typename std::multimap<TKey, TData, TCompare>::const_iterator 
			 citerator;
public:
  const TData& operator[](TKey key) const
  {
		citerator it = this->find(key);

		// if this section does not exist, throw an exception
		if (it == this->end())
		{ 
			 std::string message = std::string("Cannot find keyword '") +
						 key + "'. No such section/item exists.";
		  throw message.c_str();
		}

		// if there is more than one occurrence of 'key' no one
		// knows which should be taken. Hence throw an exception
		if (this->count(key) > 1)
			 throw "Too many identical objects.";

		return it->second;
  }

  const TData& operator()(TKey key, size_t order = 0) const
  {
		citerator it = this->lower_bound(key);

		// if the section with given order does not exist, throw an exception
		// object to prevent segmentation fault 
		if ((order > (this->count(key)-1)) || (it == this->end()))
		{
			 std::string message = std::string("Cannot find keyword '") +
						 key + "'. No such section/item exists.";
		  throw message.c_str();
		}

		while (order != 0)
		{
			 order--;
			 it++;
		}

		return (it->second);
  }

  bool present(TKey key) const
  {
		return (this->find(key) != this->end());
  }
};

//---------------------------------------------------------------------------



//---------------------------------------------------------------------------

struct ltstr
{
	 bool operator()(const char* s1, const char* s2) const 
	 { 
		  return strcmp(s1, s2) < 0; 
	 }
};

typedef my_multimap<const char*, t_store, ltstr> t_item_seq;
typedef t_item_seq* p_item_seq;

typedef my_multimap<const char *, t_item_seq, ltstr> t_section_seq;
typedef t_section_seq* p_section_seq; 

typedef std::pair<const char*, t_store> t_item;
typedef t_item* p_item;

//---------------------------------------------------------------------------

class IniHandler : public my_multimap<const char*, t_item_seq, ltstr>
{
  private:
	 std::string file_name;

/*  protected:
	 p_section_seq data;*/

  public:
	 IniHandler(const char *);
	 ~IniHandler();

	 const char *GetFileName() const
	 {
		  return file_name.c_str();
	 }

	 /*t_item_seq& operator[](const char *key)
	 {
		  return (*data)[key];
	 }

	 t_item_seq& operator()(const char *key, size_t order = 0)
	 {
		  return (*data)(key, order);
	 }

	 bool	find(const char *key)
	 {
		  return data->find(key);
	 }
	 }*/
};

//---------------------------------------------------------------------------

#endif
