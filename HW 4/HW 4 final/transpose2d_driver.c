/* test driver program for 2D block matrix transpose */
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern int transpose2d(int *a, int blocksize, MPI_Comm comm2d);

int main (int argc, char *argv[])
{
  int ROW = 0, COL = 1;  /* for readability */  
  int   numtasks, taskid;
  int i, j, n, nlocal; 
  int *alocal;
  int dims[2], periods[2];;
  MPI_Comm comm2d;
  int my2drank;
  int nrowblocks;
  int mycoords[2];
  int checkrank = 1;
  int offset;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

  if (argc != 2) {
    if (taskid == 0) {
      fprintf(stderr, "Usage: %s <n>\n", argv[0]);
      fprintf(stderr, "where n is a multiple of the square root of the number of tasks\n");
    }   
    MPI_Finalize();
    exit(0);
  }

  /* Read row/column dimension from command line */
  n = atoi(argv[1]);

  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

  dims[ROW] = dims[COL] = sqrt(numtasks);

  if (n%dims[ROW] != 0) {
    if (taskid == 0) {
      fprintf(stderr, "Usage: %s <n>\n", argv[0]);
      fprintf(stderr, "where n is a multiple of the square root of the number of tasks\n");
    }   
    MPI_Finalize();
    exit(0);
  }

  nlocal = n/dims[ROW];

  /* Create 2D Cartesian communicator */
  periods[ROW] = periods[COL] = 0;
  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &comm2d);
  MPI_Comm_rank(comm2d, &my2drank);
  MPI_Cart_coords(comm2d, my2drank, 2, mycoords);

  /*  Allocate memory for local block of A */
  alocal = (int *) malloc(nlocal*nlocal*sizeof(int));

  /*  Initialize local matrix block */
  offset = mycoords[0]*n*nlocal + mycoords[1]*nlocal + 1;
  for (i = 0; i < nlocal; i++)
      for (j = 0; j < nlocal; j++)
          alocal[i*nlocal+j] = offset + i*n + j;



  if (my2drank==checkrank) {
      printf("\noriginal matrix block for process %d\n", my2drank);
      for (i = 0; i < nlocal; i++) {
          for (j = 0; j < nlocal; j++)
              printf("%10d ", alocal[i*nlocal+j]);
          printf("\n");
      }
  }



  transpose2d(alocal, nlocal, comm2d);


 if (my2drank==checkrank) {
      printf("\nlocal transpose for process %d\n", my2drank);
      for (i = 0; i < nlocal; i++) {
          for (j = 0; j < nlocal; j++)
              printf("%10d ", alocal[i*nlocal+j]);
          printf("\n");
      }
  }



  MPI_Finalize();
  exit(0);
}
















int transpose2d(int *a, int blockdim, MPI_Comm comm2d) 
{
int *recv_array;
int temp,i,j;
MPI_Status status1,status2;


/* transposing own value*/
for (i=0;i<blockdim;i++)
for(j=0;j<i;j++)
{temp=a[j*blockdim+i];
a[j*blockdim+i]=a[i*blockdim+j];
a[i*blockdim+j]=temp;
}





/* determining own coordinates */
int mycoords[2],destcoords[2],sourcecoords[2];
int my2drank,dest_2d_rank;
int   numtasks,source_2d_rank;

MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
MPI_Comm_rank(comm2d, &my2drank);
MPI_Cart_coords(comm2d, my2drank, 2, mycoords);




/* data exchange*/
i=mycoords[0];
j=mycoords[1];

if(i!=j)
{

/*determining source and destination for my own coordinates*/ 
recv_array=(int*) malloc(blockdim*blockdim*sizeof(int));
sourcecoords[0]=j;sourcecoords[1]=i;
destcoords[0]=j; destcoords[1]=i;
MPI_Cart_rank(comm2d,sourcecoords, &source_2d_rank);
MPI_Cart_rank(comm2d,destcoords, &dest_2d_rank);



/*data exchange start */
if (i<j)
{MPI_Send(a, blockdim*blockdim, MPI_INT, dest_2d_rank, 0, comm2d);


MPI_Recv(recv_array,blockdim*blockdim, MPI_INT,source_2d_rank,1,comm2d,&status2);
}

if (i>j)
{
MPI_Recv(recv_array,blockdim*blockdim, MPI_INT,source_2d_rank,0,comm2d,&status1);

MPI_Send(a, blockdim*blockdim, MPI_INT, dest_2d_rank, 1, comm2d);
}

for (i = 0; i < blockdim; i++) 
          for (j = 0; j < blockdim; j++)
a[i*blockdim+j]=recv_array[i*blockdim+j];

free(recv_array);

}



    return(0);
}

