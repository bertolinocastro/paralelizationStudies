/****************************************************************************
 * FILE: ser_heat2D.c
 * DESCRIPTION:  
 *   Serial HEAT2D Example - C Version
 *   This example is based on a simplified 
 *   two-dimensional heat equation domain decomposition.  The initial 
 *   temperature is computed to be high in the middle of the domain and 
 *   zero at the boundaries.  The boundaries are held at zero throughout 
 *   the simulation.  During the time-stepping, an array containing two 
 *   domains is used; these domains alternate between old data and new data.
 * AUTHOR: D. Turner
 * Last Revised: 04/15/05 Blaise Barney
 ****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

struct Parms
{ 
  float cx;
  float cy;
  long int nts;
} parms = {0.1, 0.1, 10000};

int main(int argc, char *argv[])
{
	int ix, iy, iz, it;
	void inidat(), prtdat(), update();

	long unsigned tamx, tamy;

	if( argc != 3 ) exit(1);

	tamx = strtol( argv[1], NULL, 10 );
	tamy = strtol( argv[2], NULL, 10 );

	float u[2][tamx][tamy];

	// printf("Starting serial version of 2D heat example...\n");
	printf("Using [%lu][%lu] grid.\n",tamx, tamy);

	/* Initialize grid and create input file */
	// printf("Initializing grid and creating input file:");
	inidat(tamx, tamy, u);
	prtdat(tamx, tamy, u, "initial.dat");
	
	
	for (ix = 0; ix <= tamx-1; ix++) {
	   u[1][ix][0] = u[0][ix][0];
	   u[1][ix][tamy-1] = u[0][ix][tamy-1];
	   }
	for (iy = 0; iy <= tamy-1; iy++) {
	   u[1][0][iy] = u[0][0][iy];
	   u[1][tamx-1][iy] = u[0][tamx-1][iy];
	   }

	/* Iterate over all timesteps and create output file */
	// printf("Iterating over %d time steps...\n",parms.nts);





	double start = omp_get_wtime();
	iz = 0;
	for (it = 1; it <= parms.nts; it++) {
	   update(tamx, tamy, &u[iz][0][0], &u[1-iz][0][0]);
	   iz = 1 - iz;
	}
	double end = omp_get_wtime();






	printf("\n\n\n");
	printf("Tempo total: %F\n", end-start );
	// printf("Done. Created output file: ");
	prtdat(tamx, tamy, &u[iz][0][0], "final.dat");
}


/****************************************************************************
 *  subroutine update
 ****************************************************************************/
void
update(nx, ny, u1, u2)
int nx, ny;
float *u1, *u2;
{
   int ix, iy;

   for (ix = 1; ix <= nx-2; ix++) {
	  for (iy = 1; iy <= ny-2; iy++) {
		 *(u2+ix*ny+iy) = *(u1+ix*ny+iy)  + 
		 parms.cx * (*(u1+(ix+1)*ny+iy) + *(u1+(ix-1)*ny+iy) - 
		 2.0 * *(u1+ix*ny+iy)) +
		 parms.cy * (*(u1+ix*ny+iy+1) + *(u1+ix*ny+iy-1) - 
		 2.0 * *(u1+ix*ny+iy));
		 }
	  }
}

/*****************************************************************************
 *  subroutine inidat
 *****************************************************************************/
void
inidat(nx, ny, u1)
int nx, ny;
/*float u1[nx][ny];*/
float *u1;
{
   int ix, iy;

   for (ix = 0; ix <= nx-1; ix++) 
	  for (iy = 0; iy <= ny-1; iy++) 
		 *(u1+ix*ny+iy) = (float)(ix * (nx - ix - 1) * iy * (ny - iy - 1));
}

/**************************************************************************
 * subroutine prtdat
 **************************************************************************/
void
prtdat(nx, ny, u1, fnam)
int nx, ny;
float *u1;
char *fnam;
{
   int ix, iy;
   // FILE *fp;

   // fp = fopen(fnam, "w");
   // for (iy = ny-1; iy >= 0; iy--) {
   //    for (ix = 0; ix <= nx-1; ix++) {
   //      fprintf(fp, "%8.3f", *(u1+ix*ny+iy));
   //      if (ix != nx-1) {
   //         fprintf(fp, " ");
   //         }
   //     else {
   //        fprintf(fp, "\n");
   //        }
   //     }
   //  }
   // fclose(fp);
   // printf(" %s\n",fnam);

	for (iy = ny-1; iy >= 0; iy--) {
		for (ix = 0; ix <= nx-1; ix++) {
			printf("%g", *(u1+ix*ny+iy));
			if (ix != nx-1) {
				printf(" ");
			}else {
				printf("\n");
			}
		}
	}
}
