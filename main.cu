/**********************************************************************
* 
* main.cpp
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
* Description: Main body of the application (mostly the interface).
*
***********************************************************************/

//Here we define vessel3d_debug if not define
// Here defne debug if not define

#ifndef VESSEL3D_DEBUG
#define VESSEL3D_DEBUG
#endif
#ifndef DEBUG
#define DEBUG
#endif



#include <iostream>
#include "XGetopt.h"
#include<time.h>
#ifdef QT_SUPPORT
	#include <QtWidgets/QApplication>
	#include <QtWidgets/QLabel>
	#include <QtGui/QPixmap>
	#include "qtdisplay.h"
#endif

#include "ini/iniparser.h"
#include "settings.h"
#include "cpm.h"
#include "edge_detection.h"
#include<cuda.h>
#include<cuda_runtime.h>
#include "device_launch_parameters.h"
#include<device_atomic_functions.h>
#include<device_functions.h>


using namespace std;

//----------------------------------------------------------------------------------
// The function describing the use of this console application
// it is called whenever the program is called without parameters
// or with bad combination of parameters
//----------------------------------------------------------------------------------
void usage(const char *name)
{
	 cerr << "APPLICATION NAME" << endl <<
				"\t" << name << " ... 3D dynamic vessel generator  " <<
				endl << endl;
	 cerr << "SYNTAX" << endl << "\t" << name << " <options> " << endl << endl;
	cerr << "OPTIONS" << endl <<
				"\t-c <filename> ... ini file containing the description\n" <<
				"\t                  of the simulation (mandatory option)\n" <<
				"\t-p <filename> ... image with pregenerated cell population\n" <<
				"\t                  (must be compatible with the ini file)\n" <<
				"\t-h            ... this help" << endl << endl;

	 cerr << "SOME EXAMPLES:" << endl << endl;

	 cerr << "** " << "Basic generation without use of any pregenerated data:" << 
				endl <<
				"\t" << name << " -c vessel.ini" << endl << endl;
	 cerr << "** " << "Generation based on already pregenerated cell population: " <<
				endl << 
				"\t" << name << " -c config.ini -p data/cells_25.ics" <<
				endl << endl;
}

//----------------------------------------------------------------------------------
// Entry point of the program, the main() function.
//----------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
	 try {

	 // string variables with filenames
	 std::string 
				iniFilename(""),
				pregeneratedImageFilename("");

	 
	 // Read the command line input
	 if (argc < 2)
	 {
		  usage(argv[0]);
		  exit(-1);
	 }

	 // variable that controls the 'getopt' function
	 int c;

	 // read the command line options
	 while ((c = getopt(argc, argv, "c:p:h")) != -1)
	 {
			switch (c)
			{
					 case 'c': iniFilename = optarg;
								  break;
					 case 'p': pregeneratedImageFilename = optarg;
								  break;
					 case 'h': // user asks for help
					 default: usage(argv[0]);
								 exit(-1);
			}
	 }

	 if (iniFilename.empty())
	 {
		  cerr << "\nERROR: You missed to enter the configuration file! \n\n";
		  usage(argv[0]);
		  exit(-1);
	 }


	 DEBUG_REPORT("supplied configuration file: " << iniFilename.c_str());
	 IniHandler configuration(iniFilename.c_str());

	 // Initialize random seed generator (either manually or automatically)
	 unsigned int RandSeed;

	 if (configuration["cellular potts model"].present("random seed"))
	 {
		  RandSeed = (int) configuration["cellular potts model"]["random seed"];
	 }
	 else
	 {
		  RandSeed = time(NULL);
	 }

	DEBUG_REPORT("stdlib generator seeded with " << RandSeed); 
	 srand(RandSeed);

	 // Create the basic model
	 i3d::Image3d<i3d::GRAY16> img;
	 CPM model(&configuration);

	 // Set up the initial cell population
	 if (pregeneratedImageFilename.empty())
	 { 
		  model.InitializePopulation(img);
	 }
	 else
	 {
		  img= *(new i3d::Image3d<i3d::GRAY16>(pregeneratedImageFilename.c_str()));						//taking content of new object into img
		  model.ImposeInitialPopulation(img);
	 } 

	 //  Detect edges between cell population and medium
	 model.PrecomputeEdges();
/*	 EdgeDetection(model.ShowPopulation(), 
						model.GetNeighbourhood(), 
						model.GetEdgeSet());*/

	 model.Render();
	 model.StoreToFile();

	 float secreteRate = (float) configuration["pde"]["secr_rate"];
	 float decayRate = (float) configuration["pde"]["decay_rate"];
	 float diffCoeff = (float) configuration["pde"]["diff_coeff"];
	 float  diffTime = (float)configuration["pde"]["dt"];
	 float diffSpace = (float) configuration["pde"]["dx"];
	 double diffConst = (diffCoeff*diffTime) / (diffSpace*diffSpace);

	 
	 // The main application loop
#ifdef QT_SUPPORT
	 QApplication app(argc, argv);
	 MainWidget w(&model);
	 app.connect(&w, SIGNAL(SimulationCompleted(void)), SLOT(quit(void)) );
	app.exec();
#else
//----------------------------------------------------------------------------
	//Vinayak Changes
	long long  t1, t2, t3, t4, t5,t6,t7, t8;
	double secreteTime =0.0 , diffuseTime=0.0 , loopTime=0.0 , donextTime=0.0;
	 float *imgCellIDs, *cudaimgCellIDs, *cudaimgConcentration,									//Declaring Cuda Kernel Variable
		 *cudaimgAltConcentration , *imgConcentration;
	 imgCellIDs = new float[img.GetImageSize()];									
	 imgConcentration = new float[img.GetImageSize()];
	 cudaMalloc((void **)&cudaimgCellIDs, sizeof(float)*img.GetImageSize());					//Allocating Memory for Cuda Kernel Variable  
	 cudaMalloc((void **)&cudaimgAltConcentration, sizeof(float)*img.GetImageSize());			//Allocating Memory for Cuda Kernel Variable
	 cudaMalloc((void **)&cudaimgConcentration, sizeof(float)*img.GetImageSize());				//Allocating Memory for Cuda Kernel Variable
	
	 model.GetImage(imgCellIDs);														

	 cudaMemcpy(cudaimgCellIDs, imgCellIDs, sizeof(float)*img.GetImageSize(), cudaMemcpyHostToDevice);	//Copying CPU Memory to Kernel Memory(from CPU to GPU)

	 cudaMemset(cudaimgConcentration, 0.0, sizeof(float)*img.GetImageSize());							//Intializing imgConcentration to 0 in GPU 
	 cudaMemset(cudaimgAltConcentration, 0.0, sizeof(float)*img.GetImageSize());						//Intializing imgConcentration to 0 in GPU
	 int renderingPeriod = configuration["rendering"]["period"];


    t1=clock();
	 for (size_t it = 0; it < model.GetOverallDuration(); it++) 
	 {
            t2=clock();
			for (int r = 0; r < model.PdeIterations(); r++) 
			{
		        t3=clock();
				CudaSecrete << <img.GetImageSize() / 1024, 1024 >> > (cudaimgCellIDs,					//Calling Secrete Function in GPU with number of threads Per block to be 1024 as maximum capacity for this system . And passing number of blocks to be img.GetImageSize()/1024 and total number of threads are to be (img.GetImageSize()/1024)*1024. 
					cudaimgConcentration, secreteRate, diffTime, decayRate);



				t4=clock();
				CudaDiffuse << <img.GetImageSize() / 1024, 1024 >> > (cudaimgConcentration, cudaimgAltConcentration,					//Calling Diffuse function in GPU with number of threads per block to be 1 and nubmer of blocks toequal to image size distrubuted in X , Y ,Z direction according to image. 
					img.GetSizeX(), img.GetSizeY(), img.GetSizeZ(), img.GetSizeX()*img.GetSizeY(), diffConst);
				t5=clock();	
				
				secreteTime+= (double)(t4-t3)/double(CLOCKS_PER_SEC)*1000.0;
				diffuseTime+= (double)(t5-t4)/double(CLOCKS_PER_SEC)*1000.0;
			}

			//model.ShowCellVolumes();
			t6=clock();
			loopTime+=(double)(t5-t2)/double(CLOCKS_PER_SEC)*1000.0;
			
			cudaMemcpy(imgConcentration, cudaimgConcentration, sizeof(float)*img.GetImageSize(), cudaMemcpyDeviceToHost);		//Calling Cuda Memcpy to copy contents of CUda kernel memory to CPU memory .In PDe iterations only cudaimgConcentration was chaged so to copy contents of cudaimgConcentration to imgConcentration.

			model.DoNextStep(imgCellIDs, imgConcentration);

			
			cudaMemcpy(cudaimgCellIDs, imgCellIDs, sizeof(float)*img.GetImageSize(), cudaMemcpyHostToDevice);//Calling CUda Memcpy to copy contents of CPU memory to cuda Kernel Memory .In DonextStep only imgCellIDs was chaged so to copy contents of imgCellIDs to cudaimgCellIDs.
			t7=clock();
			if (it%renderingPeriod == 0) {
				model.Render(imgCellIDs);
				model.StoreToFile();

			}
            donextTime+= (double)(t7-t6)/double(CLOCKS_PER_SEC)*1000.0;
	 }
	 t8=clock();
     cout << "\n Mean GPU Time Taken By whole pRocess = " << ((t8 - t1)/double(CLOCKS_PER_SEC)*1000)/model.GetOverallDuration() << "\n";
	 cout << "\n Mean GPU TIme Taken By loop = " << loopTime / (model.GetOverallDuration()) << "\n";
	 cout << "\n Mean GPU TIme Taken By DoNextStep = " << donextTime / (model.GetOverallDuration()) << "\n";
	 cout << "\n Mean GPU TIme for secrete function = " << secretTime / (model.GetOverallDuration()*model.PdeIterations()) << "\n";
	 cout << "\n Mean GPU Time for Diffuse Funciton = " << diffuseTime / (model.GetOverallDuration()*model.PdeIterations()) << "\n";




	 delete[] imgCellIDs;
	 delete[] imgConcentration;
	 cudaFree(cudaimgCellIDs);
	 cudaFree(cudaimgAltConcentration);
	 cudaFree(cudaimgConcentration);
	 
//End of chagnes
#endif

	 } catch (std::string &e)
		  {
				cout << e << endl;
		  }
		  catch (i3d::IOException& e)
		  {
				cout << e << endl;
		  }
		  catch (i3d::InternalException& e)
		  {
				cout << e << endl;
		  }
		  catch (std::bad_alloc&)
		  {
				cout << "Not enough memory." << endl;
		  }
		  catch (...)
		  {
				cout << "System exception (2)." << endl;
		  } 

	 return 0;
}
