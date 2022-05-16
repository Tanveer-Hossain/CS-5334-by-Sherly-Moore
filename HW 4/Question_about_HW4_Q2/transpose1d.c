/* test driver program for 1D matrix transpose */
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

extern int transpose1d(int *a, int n, int blockdim, MPI_Comm comm);

int main (int argc, char *argv[])
{ 
int   numtasks, taskid;
int i, j, k ,n, nlocal; 
int *alocal;
int offset;
int checkrank = 0;

MPI_Init(&argc, &argv);
MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
MPI_Comm_rank(MPI_COMM_WORLD,&taskid);

if (argc != 2) {
   if (taskid == 0) {
     fprintf(stderr, "Usage: %s <n>\n", argv[0]);
     fprintf(stderr, "where n is a multiple of the number of tasks.\n");
   }   
   MPI_Finalize();
   exit(0);
}

/* Read row/column dimension from command line */
n = atoi(argv[1]);

if (n%numtasks) {
   if (taskid == 0) {
     fprintf(stderr, "Usage: %s <n>\n", argv[0]);
     fprintf(stderr, "where n is a multiple of the number of tasks.\n");
   }   
   MPI_Finalize();
   exit(0);
}

nlocal = n/numtasks;

/* Allocate local row block */
alocal = (int *) malloc(n*nlocal*sizeof(int));

/* Initialize local block */
offset = taskid*n*nlocal;
for (i = 0; i < nlocal; i++)
   for (j = 0; j < n; j++)
       alocal[i*n + j] = offset + i*n + j + 1;

if (taskid == checkrank) {
   printf("process %d local block before transpose\n", taskid);
   for (i = 0; i < nlocal; i++) {
       for (j = 0; j < n; j++) 
            printf("%10d ", alocal[i*n + j]);
       printf("\n");
   }
}



transpose1d(alocal, n, nlocal, MPI_COMM_WORLD);



if (taskid == checkrank) {
   printf("process %d local block after transpose\n", taskid);
   for (i = 0; i < nlocal; i++) {
       for (j = 0; j < n; j++) 
            printf("%10d ", alocal[i*n + j]);
       printf("\n");
   }
}

MPI_Finalize();
return(0);
}
















int transpose1d(int *a, int n, int blockdim, MPI_Comm comm) 

{
 
int k,i,j,temp;


for (k=0; k<n/blockdim;k++)
for (i=0;i<blockdim ; i++)
for (j=0;j<blockdim;j++)
{if(i>j)
 {temp=a[k*blockdim+i*n+j];
a[k*blockdim+i*n+j]=a[k*blockdim+j*n+i];
a[k*blockdim+j*n+i]=temp;
}
}

for (i=0;i<blockdim ; i++)
{MPI_Alltoall(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &a[i*n], blockdim, MPI_INT,comm);
MPI_Barrier( MPI_COMM_WORLD );

}



    return(0);
}
 

