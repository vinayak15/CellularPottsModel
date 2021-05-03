# CellularPottsModel
This is code for Cellula Potts model(CPM) which is implemented in C/C++. CPM is responsible for generation of cellular networks in real life based on randomness.
As CPM is very time consumin task so we have tried to optimize it by CUDA programming.
Here every time consuming function like Secrete has been implemented in CUDA  and option is given to user whether he wants to run particular function in GPU or CPU.

He can select CPU or GPU anytime in code and it becomes easier for comparison of CUDA run time with CPU runtime.

For running this code it is written in C++ and requires i3dlib and CMake for running it. 
