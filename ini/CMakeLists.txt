#######################################################################
# 
# CMakeLists.txt
#
# This file is part of VesselGen(3D)
# 
# Copyright (C) 2016 -- Centre for Biomedical Image Analysis (CBIA)
# http://cbia.fi.muni.cz/
# 
# VesselGen is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# VesselGen is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with VesselGen. If not, see <http://www.gnu.org/licenses/>.
# 
# Author: David Svoboda
# 
# Description: configuration file 
#
########################################################################

# project name
project(IniParser)

# source files
set(PARSER_SRC iniparser)
set(SOURCES 
		  ${CMAKE_CURRENT_BINARY_DIR}/lex.yy.cc 
		  ${CMAKE_CURRENT_BINARY_DIR}/${PARSER_SRC}.tab.cc 
		  iniparser.cc)

set(FLEX_SRC ${PARSER_SRC}.l)
set(BISON_SRC ${PARSER_SRC}.y)

# BISON
find_program(PROGRAM_BISON bison)
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${PARSER_SRC}.tab.cc ${CMAKE_CURRENT_BINARY_DIR}/${PARSER_SRC}.tab.hh
			  COMMAND ${PROGRAM_BISON} -d ${BISON_SRC} -o ${CMAKE_CURRENT_BINARY_DIR}/${PARSER_SRC}.tab.cc
			  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
			  DEPENDS ${BISON_SRC})

# FLEX
find_program(PROGRAM_FLEX flex)
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/lex.yy.cc
			  COMMAND ${PROGRAM_FLEX} --nounistd -o ${CMAKE_CURRENT_BINARY_DIR}/lex.yy.cc ${FLEX_SRC} 
			  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
			  DEPENDS ${FLEX_SRC} ${CMAKE_CURRENT_BINARY_DIR}/${PARSER_SRC}.tab.hh)

# final target
set(NAME "iniparser")
include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
add_library(${NAME} STATIC ${SOURCES})

install(TARGETS ${NAME} DESTINATION lib)
install(FILES iniparser.h DESTINATION include/ini)
