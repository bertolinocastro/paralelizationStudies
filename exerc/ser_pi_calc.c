/***************************************************************************
 * FILE: ser_pi_calc.c
 * DESCRIPTION:  
 *   Serial pi Calculation - C Version
 *   This program calculates pi using a "dartboard" algorithm.  See
 *   Fox et al.(1988) Solving Problems on Concurrent Processors, vol.1
 *   page 207.  
 * AUTHOR: unknown
 * REVISED: 02/23/12 Blaise Barney
***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>


void srandom (unsigned seed);  
double dboard (int darts);

#define DARTS 10000
#define ROUNDS 10000

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

        int rank, size;

        MPI_Comm_rank( MPI_COMM_WORLD, &rank );
        MPI_Comm_size( MPI_COMM_WORLD, &size );

        double pi;          
        double avepi;       
        int i;

        srandom (5);       
        avepi = 0;

        /*
         *  LOOP NAO ACEITA PARALELIZACAO VIA MPI, HA DEPENDENCIA INTERNA DA VARIAVEL avepi
         */
        for (i = rank; i < ROUNDS; i+=size) {
           pi = dboard(DARTS);
           avepi = ((avepi * i) + pi)/(i + 1);
        }

        MPI_Barrier( MPI_COMM_WORLD );

        double result[size]; result[rank] = avepi;

        //MPI_Gather( &avepi, 1, MPI_DOUBLE, &result, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );

        if(!rank){
          //  for(int i = 1; i<size; result[0]+=result[i++]);
            printf("Valor computado de pi %15.14F\n", avepi);
            printf("\nReal value of PI: 3.1415926535897 \n");
        }

    MPI_Finalize();
}


#define sqr(x)	((x)*(x))
long random(void);

double dboard(int darts)
{
   double x_coord,       
          y_coord,     
          pi,          
          r;           
   int score,          
       n;
   unsigned int cconst;
   cconst = 2 << (31 - 1);
   score = 0;

   #pragma omp paralell for schedule(dynamic) private(n,x_coord,y_coord,r) reduction(+:score)
   for (n = 1; n <= darts; n++) {
      r = (double)random()/cconst;
      x_coord = (2.0 * r) - 1.0;
      r = (double)random()/cconst;
      y_coord = (2.0 * r) - 1.0;

      if ((sqr(x_coord) + sqr(y_coord)) <= 1.0)
         score++;
    }

   pi = 4.0 * (double)score/(double)darts;
   return(pi);
} 
