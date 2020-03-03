#include "Floyd_MPI.h"

// zhou shu (22552162)
// Shikai Wang (21938451)
int main(int argc, char *argv[])
{
    struct timeval start, end;
    double delta;
    int rank,size,vexnum;// vex = the number of nodes

    // init
    PGraph PGrph;
    MPI_Comm comm;
    MPI_Init(&argc,&argv);
    gettimeofday(&start, NULL);    
    comm= MPI_COMM_WORLD;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&size);
    int offset_of_slice[size+1];
    // this array record the start offset and end offset of martix row for every threads 

    if(rank==0)
    {
        PGrph = readfile(argv[1]);

        // if size == 1, then use normal way to implement Floyd-Warshall algorithem
        
        //allocate value for different threads
        offset_of_slice[0] = 0;
        vexnum =PGrph->vexnum;
        for (int i = 1; i <= size; i++)
        {
            if (i <= PGrph->vexnum % size)
                offset_of_slice[i]=vexnum / size+offset_of_slice[i-1]+1;
            else 
                offset_of_slice[i]=vexnum / size+offset_of_slice[i-1];
        }
        
    }
    
    //split the matrix for every thread and send to every thread
    MPI_Bcast(offset_of_slice,size+1,MPI_INT,0,comm);
    MPI_Bcast(&vexnum, 1, MPI_INT, 0, comm);

    int *distance=(int*)malloc(sizeof(int)*vexnum*vexnum);// the matrix of distance
    int *distance_new=(int*)malloc(sizeof(int)*vexnum*vexnum);// the output matrix
    if (rank==0)
    {
        for (int i = 0; i < vexnum; i++)
        {
            for (int j = 0; j < vexnum; j++)
            {
                if (i==j)
                {
                    distance[i*vexnum+ j]=0;
                    continue;
                }
                
                distance[i*vexnum+j] = PGrph->matrix[i * PGrph->vexnum + j];
            }
        } 
    }
    
    // get the orginal matrix and sent to every thread
    MPI_Bcast(distance,vexnum*vexnum,MPI_INT,0,comm);    

    // Floyd-Warshall algorithm
    for (int k = 0; k < vexnum; k++)
    {
        int k_root;
        for (int index = 1; index <=size; index++)
        {
            if (k < offset_of_slice[index])
            {
                k_root = index - 1;
                break;
            }
        }
        // find which thread contain k row

        MPI_Bcast(&distance[k*vexnum], vexnum, MPI_INT, k_root, comm);

        //bcast from k_root to every thread
        for (int i = offset_of_slice[rank]; i < offset_of_slice[rank + 1]; i++)
        {
            if (distance[i*vexnum+k] != 0)
            {
                for (int j = 0; j < vexnum; j++)
                {
                    if (i!=j&&distance[k*vexnum+j] != 0 && (distance[i*vexnum+j] == 0 || distance[i*vexnum+k] + distance[k*vexnum+j] < distance[i*vexnum+j]))
                    {
                        distance[i*vexnum+j] = distance[i*vexnum+k] + distance[k*vexnum+j];
                    }
                }
            }
        }
    }

    // Gather all information
    MPI_Op MPI_Min_dis;
    MPI_Op_create((MPI_User_function *)minmum_distance, 1, &MPI_Min_dis);
    // because the disconnection is regarded as 0. Therefore, we didn't use the build-in operations
    
    MPI_Reduce(distance, distance_new, vexnum * vexnum, MPI_INT, MPI_Min_dis, 0, comm);
    gettimeofday(&end, NULL);
    delta = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
    if (rank==0)
    {
        for (int i = 0; i < vexnum; i++)
        {
            // printf("gggg");
            for (int j = 0; j < vexnum; j++)
            {
                printf("%d ", distance_new[i*vexnum+j]);
            }
            printf("\n");
        }
        outputfile(vexnum ,distance_new);
    }   
    
    

    printf("process: %d time: %12.10f\n",rank,delta);
    
    
    MPI_Finalize();
}


// this function is to reduct minmum distance
//  which was stored in every element in each thread
void minmum_distance (int *in, int*inout, int *len,MPI_Datatype *datatype)
{
    for (int i = 0; i < *len; i++)
    {
        if (inout[i]==0||(in[i]!=0 && in[i]<inout[i])) inout[i]=in[i]; 
        // in[i];
        // inout[i];
    }  
}

// test code, normal way to run Floyd

// void init(PGraph PGrph, PFloyd PFlyd)
// {
    
//     PFlyd->G = PGrph;


//     PFlyd->dist = (int *)malloc(sizeof(int) * PGrph->vexnum*PGrph->vexnum );
//     for (int i = 0; i < PGrph->vexnum; i++)
//     {
//         for (int j = 0; j < PGrph->vexnum; j++)
//         {
//             PFlyd->dist[i * PGrph->vexnum+j] = PGrph->matrix[i * PGrph->vexnum + j];
//         }
//     }   
// }

// 
// void Floyd_normal(PFloyd PFlyd)
// {
//     int n=PFlyd->G->vexnum;

//     for (int k = 0; k < PFlyd->G->vexnum; k++)
//     {
//         for (int i = 0; i < PFlyd->G->vexnum; i++)
//         {
//             if (PFlyd->dist[i*n+k] != 0)
//             {
//                 for (int j = 0; j < PFlyd->G->vexnum; j++)
//                 {
//                     if (PFlyd->dist[k*n+j] != 0 && (PFlyd->dist[i*n+j] == 0 || PFlyd->dist[i*n+k] + PFlyd->dist[k*n+j] < PFlyd->dist[i*n+j]))
//                     {
//                         PFlyd->dist[i*n+j] = PFlyd->dist[i*n+k] + PFlyd->dist[k*n+j];
//                     }
//                 }
//             }
//         }
//     }
//     for (int i = 0; i < PFlyd->G->vexnum; i++)
//     {
//         for (int j = 0; j < PFlyd->G->vexnum; j++)
//         {
//             printf(" %d ", PFlyd->dist[i*n+j]);
//         }
//         printf("\n");
//     }
//     return ;
// }
