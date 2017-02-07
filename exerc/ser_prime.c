/******************************************************************************
* FILE: ser_prime.c
* DESCRIPTION: 
*   This program generates primes. The approach taken is a "brute force"
*   method which requires increasingly greater amounts of cpu as the problem
*   size increases.  It should lend itself well to an embarassingly parallel
*   solution since each prime can be computed independently of all other
*   primes.
* AUTHOR: Blaise Barney 11/25/95 - adapted from version contributed by 
*   Richard Ng & Wong Sze Cheong during MHPCC Singapore Workshop (8/22/95).
* LAST REVISED: 04/12/05
****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>
#include <string.h>

#define exp 10

int isprime(int n) {
  int i,squareroot;
    squareroot = (int) sqrt(n);
     for (i=3; i<=squareroot; i=i+2)
        if ((n%i)==0)
           return 0;
     return 1;
}


int main(int argc, char *argv[]){
  int pc, LIMIT = ( argc == 2 )? strtol(argv[1],NULL,10) : pow( 2, exp );

  MPI_Init(&argc, &argv);
    int rank, size;

    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );

    pc=(!rank)?4:0;

    int passo = (LIMIT-11)/size;
    int excesso = (LIMIT-11)%size;

    int inicio = 11+passo*rank;
    int fim=(rank==size-1)?inicio+passo+excesso:inicio+passo;

    int pass = size*2;

    double start = omp_get_wtime();

    omp_set_num_threads( omp_get_num_procs() );
    for (int j = 0; j < 20; ++j)
    {
    // #pragma omp parallel for schedule( dynamic ) reduction(+:pc) /* Best */
    //   for( int i=(!(inicio%2))?inicio+1:inicio; i<fim; i+=2)
    //       if(isprime(i))
    //           pc++;


    omp_set_num_threads( omp_get_num_procs() );
    #pragma omp parallel for schedule(dynamic) reduction(+:pc)/* Worse */
    for (int i = 11+rank*2; i < LIMIT; i+=pass)
      if(isprime(i))
        pc++;

    }
    

    

    
    double end = omp_get_wtime();
    MPI_Barrier( MPI_COMM_WORLD );

    int result[size];
    double tempo[size], mytmp=end-start;

    MPI_Gather( &pc, 1, MPI_INT, &result, 1, MPI_INT, 0, MPI_COMM_WORLD );
    MPI_Gather( &mytmp, 1, MPI_DOUBLE, &tempo, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );

    if( !rank ){
      double temp=mytmp;
      for (int i = 1; i < size; ++i)
        pc += result[i],
        temp+= tempo[i];

      printf("Total: %d tempo %F\n", pc, temp/size);
    }

  MPI_Finalize();

  return 0;
} 
