/**********************************************************************
* 
* iniparser.y
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

%{
#include <stdio.h>
#include <iostream>
#include "iniparser.h"

#define YYERROR_VERBOSE

int yyerror(p_section_seq *data, const char *str);
int yylex();

extern int yylineno;
%}

%union 
{
	int					int_type;
	char					*str_type;
	float					float_type;
	t_couple 			couple_type;
	t_triplet 			triplet_type;
	t_store				store_type;
	p_item				item_ptype;
	p_item_seq			item_seq_ptype;
	p_section_seq		section_seq_ptype;
}

/* tokens */
%token UNKNOWN_TOKEN
%token <str_type>	   NAMES	
%token <int_type>		INUMBER
%token <float_type>	FNUMBER

/* nonterminal types */
%type <float_type>				number
%type <triplet_type>				triplet
%type <couple_type>				couple
%type <store_type>				keyvalue
%type <item_ptype>				item	
%type <item_seq_ptype>			item_seq
%type <section_seq_ptype>		section_seq

%parse-param {p_section_seq *data}

%%

/* grammar rules */
root: section_seq 
	 { 
	 memcpy(data, &$1, sizeof(void*));
	 }
	 | '\n' section_seq 
	 { 
	 memcpy(data, &$2, sizeof(void*));
	 }
;

section_seq:	'[' NAMES ']' '\n'
	 { 
	 	#ifdef DEBUG
			std::cout << "[" << $2 << "]\t(section name)" << std::endl;
		#endif
	 } item_seq
	 {
		$$ = new t_section_seq;
		$$->insert(std::pair<const char *, t_item_seq>($2, *$6));
		delete $6;
	 }
	 		| section_seq '[' NAMES ']' '\n' 
			{
				#ifdef DEBUG
					std::cout << "[" << $3 << "]\t(section name)" << std::endl;
				#endif
			} item_seq
			{
				$$ = $1;
		 		$$->insert(std::pair<const char *, t_item_seq>($3, *$7));
				delete $7;
			}
;

item_seq:	item
	 { 
		 $$ = new t_item_seq; 
		 
		 if ($1 != NULL)
		 { 
		 	//(*$$)[$1->first] = $1->second; 
			$$->insert(*$1);
			delete $1; 
		 }
	 }
	 | item_seq item 
	 { 
		$$ = $1;
		if  ($2 != NULL)
		{
			//(*$$)[$2->first] = $2->second; 
			$$->insert(*$2);
			delete $2; 
		}
	 }
;

item: '\n' /* empty line */ 
	 { 
	 	$$ = NULL; 
	 }
	 | NAMES 
	 {
	 	#ifdef DEBUG
	   	std::cout << "'" << $1 << "' = "; 
	 	#endif
	 }
	 '=' keyvalue '\n'
	 { 
	 	$$ = new t_item($1,$4); 
	 }
;

keyvalue: 	triplet
		  		{ 
					$$ = $1; 
					#ifdef DEBUG 
						std::cout << $1 << "\t(triplet)" << std::endl; 
					#endif
				}
		  		| couple 
		  		{ 
					$$ = $1; 
					#ifdef DEBUG 
						std::cout << $1 << "\t(couple)" << std::endl; 
					#endif
				}
		  		| INUMBER
		  		{ 
					$$ = $1; 
					#ifdef DEBUG 
						std::cout << $1 << "\t(integer number)" << std::endl; 
					#endif
				}
		  		| FNUMBER
		  		{ 
					$$ = $1; 
					#ifdef DEBUG 
						std::cout << $1 << "\t(floating point number)" << std::endl; 
					#endif
				}
				| NAMES 
		  		{ 
					$$ = $1; 
					#ifdef DEBUG 
						std::cout << "'" << $1 << "'\t(string)" << std::endl; 
					#endif
				}
				| NAMES '(' NAMES ')'
		  		{ 
					char tmp[MAX_STR_LEN];
					strcpy(tmp, $1);
					strcat(tmp, "(");
					strcat(tmp, $3);
					strcat(tmp, ")");

					$$ = tmp; 
					#ifdef DEBUG 
						std::cout << "'" << $1 << " (" << $3 << ")" << 
										"'\t(string)" << std::endl; 
					#endif
				}
;

triplet: '(' number ',' number ',' number ')'
		{ 
			$$.x = $2;
			$$.y = $4;
			$$.z = $6;
		}
;

couple: '(' number ',' number ')' 
		{
			$$.x = $2;
			$$.y = $4;
		}
;

number: 	INUMBER 
		{ 
			$$ = $1; 
		}
		| FNUMBER 
		{ 
			$$ = $1; 
		}
;

/* end of grammar rules */
%%

int yyerror(p_section_seq *data, const char *str)
{
	printf ("INI file parse error: line %d: \n -> %s\n",yylineno,str);
	return 0;
}

extern "C" int yywrap(void)
{
	return 1;
}
