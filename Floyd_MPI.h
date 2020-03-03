#ifndef _FLOYD_MPI_
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include<mpi.h>
# include<sys/time.h>


#define MAX 4096

typedef struct 
{
    int vexnum;    
    int *matrix;
}Graph, *PGraph;

PGraph readfile(char*);
void print_gragh(PGraph);
void Floyd_MPI(int **,int,int,int*,int,MPI_Comm);
void minmum_distance (int *in, int*inout, int*len,MPI_Datatype *datatype);
void outputfile (int,int*);

// typedef struct
// {
//     PGraph G;
//     int *dist;
// }Floyd,*PFloyd; // used in normal Floyd algorithm

// void init(PGraph, PFloyd);
// void Floyd_normal(PFloyd);

#endif
