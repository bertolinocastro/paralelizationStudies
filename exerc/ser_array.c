/******************************************************************************
* FILE: ser_array.c
* DESCRIPTION: 
*   Serial Example - Array Assignment - C Version
*   In this simple example, an array is initialized and values assigned.
* AUTHOR: Blaise Barney
* LAST REVISED:  04/15/05
****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>

#define	ARRAYSIZE	1000000

int main(int argc, char *argv[]){

	MPI_Init(&argc, &argv);
	
		int rank, size, tot_proc;
		MPI_Comm_rank( MPI_COMM_WORLD, &rank ),
		MPI_Comm_size( MPI_COMM_WORLD, &size );

		int     i; 			/* loop variable */
		float	data[ARRAYSIZE]; 	/* the intial array */
		double start;

		//if( !rank )
			start = MPI_Wtime();
		if( !rank )
			tot_proc = omp_get_num_procs();

		#pragma omp parallel num_threads( tot_proc )
		{	
			printf("Sou o thread nº %d do processo %d\n", omp_get_thread_num(), rank);

			#pragma omp for private(i) schedule( guided )
				for(i=0; i<ARRAYSIZE; i++)
				  data[i] =  i * 1.0,
				  data[i] = data[i] + i * 1.0;
		}

		//if( !rank )
			printf("Tempo decorrido %F\n", MPI_Wtime()-start);

	MPI_Finalize();

	return 0;
}
