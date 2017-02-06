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

#define exp 5

int isprime(int n) {
  int i,squareroot;
  if (n>10) {
     squareroot = (int) sqrt(n);
     for (i=3; i<=squareroot; i=i+2)
        if ((n%i)==0)
           return 0;
     return 1;
     }
  else
   return 0;
}


int main(int argc, char *argv[]){               
  int pc, LIMIT = pow( 2, exp );

  pc=4;

  MPI_Init(&argc, &argv);
    int rank, size;

    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );

    int passo = (LIMIT-11)/size;
    
    MPI_Status status[size];
    MPI_Request request[size];
  
  //if( !rank ){
      omp_set_num_threads( omp_get_num_procs() );
      #pragma omp parallel for schedule( dynamic ) shared(pc)
        for(int i=11 + passo*rank; i<11+passo*(rank+1); i+=2)
          {
            printf( "Sou o proce %d, thread %d de passo %i\n", rank, omp_get_thread_num(), i );
            if(isprime(i))
              #pragma omp atomic
                pc++;
          }
      
  //}
      MPI_Barrier( MPI_COMM_WORLD );

      int result[size];

      MPI_Gather( &pc, 1, MPI_INT, &result, 1, MPI_INT, 0, MPI_COMM_WORLD );

      if( !rank ){
        for (int i = 1; i < size; ++i)
          pc += result[i];

        printf("Total: %d\n", pc);
      }
      


  MPI_Finalize();

  return 0;
} 
