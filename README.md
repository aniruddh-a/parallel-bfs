## Parallel Breadth First Search On Shared Memory                                                                   
++++++++++++++++++++++++++++++++++                                                           
By Aniruddh Ramesh Adkar and Anuj Sharma

All source code used in the project are for curriculum use only.                              
++++++++++++++++++++++++++++++++++

Project : Parallel breadth first search on shared memory system.                                

Requirements:                                                                                  

1) Any operating system supporting gcc compiler and installed openmpi environment along with boost library
2) We have tested with gcc 4.8.4 and mpic++ 4.8.4, boost 1.5.9
3) System with multiple cores or any HPC supporting clustered nodes for parallel computing.     


Compile : 
mpic++ -std=c++11 -O3 parallel_bfs.cpp -o parallel_bfs.o

Execute :
mpiexec -n <no-of-processors> parallel_bfs.o  1000 100 0 

Usage : 
Usage: <program_name> no-of-vertices no-of-vertices-per-graph source-vertex 

Change logs:                                                                                     

WED, Oct 30, 2015                                                                              
Git repo initialized                                                                           
Project branched                                                                                  

*****************

MON, Nov 09, 2015                                                                               
Implemented sequential form of breath first search.                                            
Parallel algorithm using 1D partitioning in progress and sequential form will help in comparing results and to gain insight in computation speedup achieved by parallelism. Initially generated data for 60,000 graph vertices.       

added sequential_bfs.c                                                                         
added graphData.txt                                                                           

==================================                                                                 
