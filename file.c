#include "Floyd_MPI.h"

PGraph readfile(char* filename)
{
    FILE *fp;
    PGraph G;
    G=(PGraph)malloc(sizeof(Graph));

    if ((fp = fopen(filename, "rb")) == NULL)//open file
    {
        printf("cannot open file\n");
        exit(0);
    }
    fread(&G->vexnum,sizeof(int),1,fp);
    G->matrix = (int*)malloc(sizeof(int) * G->vexnum*G->vexnum);
    fread(G->matrix, sizeof(int), G->vexnum * G->vexnum, fp);
    return(G);
}

void print_gragh(PGraph G)
{
    printf("%d\n",G->vexnum);

    for (size_t i = 0; i < G->vexnum; i++)
    {
        for (size_t j = 0; j < G->vexnum; j++)
        {
            printf("%d ",G->matrix[i*G->vexnum+j]);
        }
        printf("\n");
    }
}
void outputfile(int num_node,int * distance)
{

    FILE *fw;
    FILE *fw_b;
    
    
    char txt_name[16];
    char bin_name[16];
    sprintf(txt_name, "%d.out.txt", num_node);
    sprintf(bin_name, "%d.out.bin", num_node);

// binary file
    fw_b = fopen(bin_name, "wb+");
    fw = fopen(txt_name, "w");
    if(fw == NULL||fw_b==NULL)
    {
        printf("cannot write in file");
    }
    fwrite(&num_node,sizeof(int),1,fw_b);
    fwrite(distance,sizeof(int),num_node*num_node,fw_b);
    fclose(fw_b);

// txt file

    
    fprintf(fw,"%d\n",num_node);
    for (size_t i = 0; i < num_node; i++)
    {
        
        for (size_t j = 0; j < num_node; j++)
        {
            fprintf(fw,"%d ",distance[i*num_node+j]);
        }
        
        fputc('\n',fw);
        
    }
    
    fclose(fw);
}
