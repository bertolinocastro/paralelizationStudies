/***************************************************************************
 * FILE: ser_wave.c
 * DESCRIPTION:
 *   Serial Concurrent Wave Equation - C Version
 *   This program implements the concurrent wave equation described 
 *   in Chapter 5 of Fox et al., 1988, Solving Problems on Concurrent
 *   Processors, vol 1. 
 * AUTHOR: unknown
 * LAST REVISED:  04/15/05 Blaise Barney
***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <mpi.h>
#include <string.h>

#define MAXPOINTS 100
#define MAXSTEPS  1000000
#define MINPOINTS 20
#define PI 3.14159265

#define TPOINTS   100
#define NSTEPS    1000000


void init_param(void);
void init_line(void);
void update (void);
void printfinal (void);

const int nsteps = NSTEPS,      /* number of time steps */
    tpoints = TPOINTS,          /* total points along string */
    rcode;                  	/* generic return code */
double values[MAXPOINTS+2], 	/* values at time t */
       oldval[MAXPOINTS+2], 	/* values at time (t-dt) */
       newval[MAXPOINTS+2]; 	/* values at time (t+dt) */

size_t mmcpytam;

const double    dtime = 0.3,
                c     = 1.0,
                dx    = 1.0;

double tau, sqtau;

/***************************************************************************
 *	Obtains input values from user
 ***************************************************************************/
void init_param(void)
   {
   // char tchar[8];

   // /* set number of points, number of iterations */
   // tpoints = 0;
   // nsteps = 0;
   // while ((tpoints < MINPOINTS) || (tpoints > MAXPOINTS)) {
   //    printf("Enter number of points along vibrating string [%d-%d]: ",
   //           MINPOINTS, MAXPOINTS);
   //    scanf("%s", tchar);
   //    tpoints = atoi(tchar);
   //    if ((tpoints < MINPOINTS) || (tpoints > MAXPOINTS))
   //       printf("Invalid. Please enter value between %d and %d\n", 
   //               MINPOINTS, MAXPOINTS);
   //    }
   // while ((nsteps < 1) || (nsteps > MAXSTEPS)) {
   //    printf("Enter number of time steps [1-%d]: ", MAXSTEPS);
   //    scanf("%s", tchar);
   //    nsteps = atoi(tchar);
   //    if ((nsteps < 1) || (nsteps > MAXSTEPS))
   //       printf("Invalid. Please enter value between 1 and %d\n", MAXSTEPS);
   //    }

   // printf("Using points = %d, steps = %d\n", tpoints, nsteps);

   

   }

/***************************************************************************
 *     Initialize points on line
 **************************************************************************/
void init_line(void){
    int j; double x, k = 0.0;
    /*register const */double fac = 2.0 * PI,
                          tmp = tpoints - 1;

    //#pragma omp parallel for schedule(static,1) private(j) ordered
    for (j = 1; j <= tpoints; j++) {
        //#pragma omp ordered
        x = k/tmp;
        values[j] = sin (fac * x);
        //#pragma omp ordered
        k = k + 1.0;
    }

    //#pragma omp parallel for private(i) schedule(dynamic)
    // for (i = 1; i <= tpoints; i++)
    //     oldval[i] = values[i];

    // double tes1, tes2;
    // tes1 = omp_get_wtime();
    memcpy(
        &oldval[1],
        &values[1],
        mmcpytam
    );
    // tes2 = omp_get_wtime();
    // printf("Custo copia vetores %f\n", tes2-tes1);
    // exit(1);
}
//void *memcpy(void *str1, const void *str2, size_t n)
/***************************************************************************
 *      Calculate new values using wave equation
 **************************************************************************/

void update(){
    int i, j;
    //#pragma omp parallel for schedule(dynamic) private(i,j)
    for (i = 1; i<= nsteps; i++) {
        newval[1] = 0.0;
        //#pragma omp parallel for schedule(static,1) private(j)
        for (j = 2; j < tpoints; j+=5) {
            // if( j == 1 || j == tpoints )
            //     newval[j] = 0.0;
            // else
                
                newval[j] = (2.0 * values[j]) - oldval[j] 
               + (sqtau * (values[j-1] - (2.0 * values[j]) + values[j+1]));

                newval[j+1] = (2.0 * values[j+1]) - oldval[j+1] 
               + (sqtau * (values[j] - (2.0 * values[j+1]) + values[j+2]));

                newval[j+2] = (2.0 * values[j+2]) - oldval[j+2] 
               + (sqtau * (values[j+1] - (2.0 * values[j+2]) + values[j+3]));

                newval[j+3] = (2.0 * values[j+3]) - oldval[j+3] 
               + (sqtau * (values[j+2] - (2.0 * values[j+3]) + values[j+4]));


                newval[j+4] = (2.0 * values[j+4]) - oldval[j+4] 
               + (sqtau * (values[j+3] - (2.0 * values[j+4]) + values[j+5]));
        }
        newval[tpoints] = 0.0;

        // for (j = 1; j <= tpoints; j++) {
        //     oldval[j] = values[j];
        //     values[j] = newval[j];
        // }


        memcpy(
            &oldval[1],
            &values[1],
            mmcpytam
        );
        memcpy(
            &values[1],
            &newval[1],
            mmcpytam
        );
    }
    /*for (int ij = 1; ij<= nsteps*tpoints; ij++) {
        int ijj = (ij%tpoints);
        if ((ijj == 1) || (ijj  == tpoints))
                newval[ijj] = 0.0;
            else
                do_math(ijj);
    }*/
}

void printfinal()
   {
   int i;
   for (i = 1; i <= tpoints; i++) {
      printf("%6.4f ", values[i]);
      if (i%10 == 0)
         printf("\n");
      }
   }

int main(int argc, char *argv[])
{
/*
int left, right;
*/
    mmcpytam = tpoints*sizeof(double);

    tau = (c * dtime / dx);
    sqtau = tau * tau;

    omp_set_num_threads( omp_get_max_threads() );

      init_param();
      
      
    double start = omp_get_wtime();
    for(int i = 0; i < 1000; ++i)
        init_line(),
        update();
    double end = omp_get_wtime();

      

      printfinal();

      printf("\n\n\nTempo gasto %.7F\n", end-start);

   return 0;
}
