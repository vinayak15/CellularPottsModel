# CellularPottsModel
This is code for Cellula Potts model(CPM) which is implemented in C/C++. CPM is responsible for generation of cellular networks in real life based on randomness.
As CPM is very time consumin task so we have tried to optimize it by CUDA programming.
Here every time consuming function like Secrete has been implemented in CUDA  and option is given to user whether he wants to run particular function in GPU or CPU.

He can select CPU or GPU anytime in code and it becomes easier for comparison of CUDA run time with CPU runtime.

For running this code it is written in C++ and requires i3dlib and CMake for running it. 


Function which are important for understanding of code
1) i3dlib is Center of Biomedical for Image Analysis lab library which can be downloaded from https://cbia.fi.muni.cz/software/i3d-library.html

2) Folder Data - It consists of input image of varying sizes.
3) main.cu - Program starts at this file , so we create a CPM model , get input image, output file  and all other functions are called and handled here.
Some time intensive functions which were handled with great efficiency.
Secrete, Diffuse, Do next step.

We observed an improvement of 10% improvement in timing when the whole program ran in CUDA as compared to CPU.
Here different functions can be tested in run time on CPU and GPU , and cases are handled if some parts of the program run in CUDA and other parts in CPU.
This is done to make sure we understand whether parallelizing function has improved timing or not.
