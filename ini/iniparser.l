/**********************************************************************
* 
* iniparser.l
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

%top{
#include <stdlib.h>
#include "iniparser.h"
#include "iniparser.tab.hh"

#define YY_NO_UNISTD_H
}

DIGIT		[0-9]
NAMES		([^[:space:]=\[\]\(\),;]|[^[:space:]=\[\]\(\),;][^=\[\]\(\),;\n]*[^[:space:]=\[\]\(\),;])
/* NAMES = string including space expect at the begining and at the end. It must not include characters like these six: ][)(,;
*/

%option yylineno
%option nounput
%option never-interactive

%%

^";".*"\n"		/* eat up comments including the new line symbol */
[+-]?{DIGIT}+	yylval.int_type = atoi(yytext); return INUMBER;
[+-]?{DIGIT}+"."{DIGIT}*	yylval.float_type = (float) atof(yytext); return FNUMBER;
{NAMES}	yylval.str_type = strdup(yytext); return NAMES;
[ \t]+		/* eliminate whitespaces */
[ \t\n]*"\n"[ \t\n]*	return '\n'; /* reduce whitespaces and new line
												 to just one new line */
"["			return '[';
"]"			return ']';
"("			return '(';
")"			return ')';
"="			return '=';
","			return ',';
.				yylval.str_type = strdup(yytext); return UNKNOWN_TOKEN;

%%

