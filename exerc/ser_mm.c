/******************************************************************************
* FILE: ser_mm.c
* DESCRIPTION:  
*   Serial Matrix Multiply - C Version
* AUTHOR: Blaise Barney
* LAST REVISED: 04/12/05
******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define NRA 62 			/* number of rows in matrix A */
#define NCA 15			/* number of columns in matrix A */
#define NCB 7   		/* number of columns in matrix B */

int main(int argc, char *argv[]){

    int    i, j, k;			/* misc */
    double a[NRA][NCA], 		/* matrix A to be multiplied */
           b[NCA][NCB],      	/* matrix B to be multiplied */
           c[NRA][NCB];		/* result matrix C */

    double start = omp_get_wtime();
    omp_set_num_threads( omp_get_max_threads() );
    #pragma omp parallel
    {
        #pragma omp for collapse(2) private(i,j) schedule(dynamic)
            for (i=0; i<NRA; i++)
               for (j=0; j<NCA; j++)
                  a[i][j]= i+j;

        #pragma omp for collapse(2) private(i,j) schedule(dynamic)
            for (i=0; i<NCA; i++)
               for (j=0; j<NCB; j++)
                  b[i][j]= i*j;

        #pragma omp for collapse(2) private(i,j) schedule(dynamic)
            for(i=0;i<NRA;i++)
                for(j=0;j<NCB;j++)
                  c[i][j] = 0.0;
    
        
        #pragma omp for collapse(3) private(i,j,k) schedule(dynamic)
        for(i=0;i<NRA;i++)
           for(j=0;j<NCB;j++)
              for(k=0;k<NCA;k++)
                 c[i][j]+= a[i][k] * b[k][j];

    }
    printf("Tempo: %F\n", omp_get_wtime()-start);
    return 0;
}
