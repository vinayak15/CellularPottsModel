/**********************************************************************
* 
* rnd_generators.cpp
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
* Author: Vladimir Ulman
* 
* Description: Various random generators (using GSL software package).
*
***********************************************************************/

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <float.h>

#include "rnd_generators.h"
#include "settings.h"

/*
 * Based on Pierrre L'Ecuyer, http://www.iro.umontreal.ca/~lecuyer/,
 * well performing random number generators also available in the GSL are:
 *
 * gsl_rng_mt19937
 * gsl_rng_taus2
 *
 * In case, one would like to give them a try, change lines
 *
 *		randState = gsl_rng_alloc(gsl_rng_default);
 * to, for example,
 *		randState = gsl_rng_alloc(gsl_rng_mt19937);
 *
 * Both "super generator", probably, need no extra seeding.
 * They should do seed themselves somehow...
 * http://www.gnu.org/software/gsl/manual/html_node/Random-number-generator-algorithms.html
 */

float GetRandomGauss(const float mean, const float sigma) {
	static gsl_rng *randState=NULL;

	//do we need to init random generator?
	if (randState == NULL) {
		//yes, we do:
		//create instance of the generator and seed it
		randState = gsl_rng_alloc(gsl_rng_default);
		unsigned long s=-1 * (int) time(NULL);
		DEBUG_REPORT("randomness started with seed " << s);
		gsl_rng_set(randState,s);
	}

	return ( gsl_ran_gaussian(randState, sigma) + mean );
}


float GetRandomUniform(const float A, const float B) {
	static gsl_rng *randState=NULL;

	//do we need to init random generator?
	if (randState == NULL) {
		//yes, we do:
		//create instance of the generator and seed it
		randState = gsl_rng_alloc(gsl_rng_default);
		unsigned long s=-1 * (int) time(NULL);
		DEBUG_REPORT("randomness started with seed " << s);
		gsl_rng_set(randState,s);
	}

	return ( gsl_ran_flat(randState, A,B) );
}

size_t GetRandomUniform(const size_t A, const size_t B) 
{
	 return (size_t) floorf(GetRandomUniform((const float) A, (const float) B+(1.0f-0.00001f)));
}

unsigned int GetRandomPoisson(const float mean) {
	static gsl_rng *randState=NULL;

	//do we need to init random generator?
	if (randState == NULL) {
		//yes, we do:
		//create instance of the generator and seed it
		randState = gsl_rng_alloc(gsl_rng_default);
		unsigned long s=-1 * (int) time(NULL);
		DEBUG_REPORT("randomness started with seed " << s);
		gsl_rng_set(randState,s);
	}

	return ( gsl_ran_poisson(randState, mean) );
}


template <typename MT>
void SuggestBrownianVector(i3d::Vector3d<MT> &v, const MT step)
{
	 // wrapper for low level function
	 SuggestBrownianVector(v, i3d::Vector3d<MT>(step, step, step));
}

template <typename MT>
void SuggestBrownianVector(i3d::Vector3d<MT> &v, const i3d::Vector3d<MT> step)
{
	//1.6f is a correction coefficient such that the mean size
	//of the generated vectors is 'step', variance of generated
	//vectors will cca 0.42 the mean - that's just a fact :-)
	/*v.x=static_cast<MT>((float)step/1.6f * GetRandomGauss(0,1));
	v.y=static_cast<MT>((float)step/1.6f * GetRandomGauss(0,1));
	v.z=static_cast<MT>((float)step/1.6f * GetRandomGauss(0,1));*/

	// for 2x bigger cells, we will force the Brownian vector to be stronger
	v.x=static_cast<MT>(1.0f * (float)step.x/1.6f * GetRandomGauss(0,1));
	v.y=static_cast<MT>(1.0f * (float)step.y/1.6f * GetRandomGauss(0,1));
	v.z=static_cast<MT>(1.0f * (float)step.z/1.6f * GetRandomGauss(0,1));
}

template void SuggestBrownianVector(i3d::Vector3d<float> &v, const float step);
template void SuggestBrownianVector(i3d::Vector3d<float> &v, const i3d::Vector3d<float> step);


float GetRandomRotationMultimodal(void)
{
	/*
	 * The probability function is:
	 * plot [-35:35] (35-abs(x))/35*0.6  *0.5*cos(x*0.571198)+ (35-abs(x))/35  *(0.5+(1-0.6)/2);
	 * 
	 * where 35 is the range in one direction,
	 * 0.571198 converts deg to rad such that 11deg is 2PI,
	 * 0.6 shrinks the curve such that it does not touch zero,
	 * (35-abs(x)/35 forces the curve to approach zero at tails
	 */

	//some tunning macros
	#define MAXDEG		35
	#define PERIODDEG	11.f
	#define SPREAD  	0.6f

	//derived helper values
	#define PERIODRAD	(6.28318f/PERIODDEG)
	#define SHIFT		((1.f-SPREAD)/2.f)

	//the function helper look up table,
	//basically, it is a distribution function
	static float LUT[2*MAXDEG +1]={-5.f};

	if (LUT[0] == -5.f) {
		LUT[0]=0.f;
		//LUT is not initialized, do it now
		for (int x=-MAXDEG+1; x <= MAXDEG; ++x) {
			//somewhat a probability of value x
			float R=((float)MAXDEG-std::fabs(x))/(float)MAXDEG;
			R*=SPREAD*0.5f*cos((float)x*PERIODRAD) + 0.5f + SHIFT;

			LUT[x+MAXDEG]=R+LUT[x+MAXDEG-1];
			//DEBUG_REPORT(x << " " << LUT[x+MAXDEG]);
		}

		DEBUG_REPORT("LUT initialized, sum=" << LUT[2*MAXDEG]);
	}

	//proceed using the "inversion technique" for generating random variates
	float val=GetRandomUniform(0,LUT[2*MAXDEG]);

	int ang=-MAXDEG;
	while ((ang <= MAXDEG) && (val > LUT[ang+MAXDEG])) ++ang;

	//consistency check
	if (ang > MAXDEG) {
		ang=MAXDEG;
		DEBUG_REPORT("WTF? LUT table was overwritten or something...");
	}

	//convert ang to float and radians
	return((float)ang *0.0174532f);
}

