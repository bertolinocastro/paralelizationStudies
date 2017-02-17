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
#include <mpi.h>

struct Parms{
	float cx;
	float cy;
	long int nts;
} parms = {0.1, 0.1, 10000};

int main(int argc, char *argv[]){
	
	float *u;
	long unsigned ix, iy, it;
	int  z;
	void inidat(), prtdat();
	long unsigned tamx, tamy;

	if( argc != 3 ) exit(1);

	tamx = strtol( argv[1], NULL, 10 );
	tamy = strtol( argv[2], NULL, 10 );

	// tamx*tamy

	u = malloc( 2*tamx*tamy*sizeof(*u) );
	if(!u) exit(1);
	// for (int i = 0; i < 2; ++i){
	// 	u[i] = malloc( tamx*szeof(**u) );
	// 	if(!u[i]) exit(1);
	// 	for (int j = 0; j < tamx; ++j){
	// 		u[i][j] = malloc( tamy*szeof(***u) );
	// 		if(!u[i][j]) exit(1);
	// 	}
	// }












	printf("Using [%lu][%lu] grid.\n",tamx, tamy);

	inidat(tamx, tamy, u);
	 prtdat(tamx, tamy, u);
 // prtdat(tamx, tamy, (u+tamx*tamy));

	for(ix = 0; ix < tamy; ix++){
		u[(tamx*tamy)+ix*tamx] = u[ix*tamx];
		u[(tamx*tamy)+ix*tamx+tamx-1] = u[ix*tamx+tamx-1];
	}

	for(iy = 0; iy < tamx; iy++){
		u[(tamx*tamy)+iy] = u[iy];
		u[(tamx*tamy)+(tamy-1)*tamx+iy] = u[(tamy-1)*tamx+iy];
	}
 
 // prtdat(tamx, tamy, (u+tamx*tamy));


// tamx colunas 
// tamy linhas

// tamx*tamy celulas // por matrz

// tamx*i // escolhe linha

// tamy*j // escolhe coluna









	double start = omp_get_wtime();
	z = 0;
	long unsigned i, j;
	for(it = 1; it <= parms.nts; it++){
		#pragma omp parallel for num_threads(3) firstprivate(z) private(i,j) schedule(dynamic)
		//{
		// 	int th = omp_get_thread_num();
		// 	int tt = omp_get_num_threads();
			// printf("Th %d TT %d\n", th, tt);
				 // #pragma omp for schedule(static) private(ix,iy, ij)
					// for (ix = 1; ix < tamx-1; ix++) {
					//     for (iy = 1; iy < tamy-1; iy++) {
					//         u[!z][ix][iy] = u[z][ix][iy]  + 
					//         parms.cx * (u[z][ix+1][iy] + u[z][ix-1][iy] - 
					//         2.0 * u[z][ix][iy]) +
					//         parms.cy * (u[z][ix][iy+1] + u[z][ix][iy-1] - 
					//         2.0 * u[z][ix][iy]);
					//     }
					// }




			// for(i ; i < x;++i ){
			// 	for (j; j < y; ++j){
			// 		(m + z*(tamx*tamy) + i*tamx + j) posicao geral do espaco
			// 	}
			// }


			for(i = 1; i < tamx-1; ++i){
				for (j = 1; j < tamy-1; ++j){
					*(u + (!z)*tamx*tamy + i*tamx + j) = *(u + z*tamx*tamy + i*tamx + j) + 
						parms.cx * (*(u + z*tamx*tamy + (i)*tamx + j+1) + *(u + z*tamx*tamy + (i)*tamx + j+1) - 2.0*(*(u + z*tamx*tamy + i*tamx + j))) + 
						parms.cy * (*(u + z*tamx*tamy + (i+1)*tamx + j) + *(u + z*tamx*tamy + (i-1)*tamx + j) - 2.0*(*(u + z*tamx*tamy + i*tamx + j)));
				}
			}






					// for(ij = th; ij < xy; ij+=tt){
					// 	ix = ij/limx+1;
					// 	iy = ij%limy+1;
					// 	//printf("th %3d ij: %5lu ix: %5lu iy: %5lu\n", th, ij, ix, iy);
					// 	u[!z*(tamx*tamy)+ix*tamy+iy] = u[z*(tamx*tamy)+ix*tamx+iy]  + 
					// 	    parms.cx * (u[z*(tamx*tamy)+(ix+1)*tamx+iy] + u[z*(tamx*tamy)+(ix-1)*tamx+iy] - 
					// 	    2.0 * u[z*(tamx*tamy)+ix*tamx+iy]) +
					// 	    parms.cy * (u[z*(tamx*tamy)+ix*tamx+iy+1] + u[z*(tamx*tamy)+ix*tamx+iy-1] - 
					// 	    2.0 * u[z*(tamx*tamy)+ix*tamx+iy]);
			// }
		//}	
		z = !z;
	}
	
	
	double end = omp_get_wtime();












	printf("\n\n\n");
	prtdat(tamx, tamy, (u + z*(tamx*tamy)));

	printf("Tempo total: %F\n", end-start );

	// for (int i = 0; i < 2; ++i){
	// 	for (int j = 0; j < tamx; ++j)
	// 		free( u[i][j] );
	// 	free( u[i] );
	// }
	free( u );
	return 0;
}












void inidat( tamx, tamy, u1 )
long unsigned tamx, tamy;
float *u1;
{
	int ix, iy;

	for (ix = 0; ix < tamx; ix++) 
		for (iy = 0; iy < tamy; iy++) 
			u1[ix*tamx+iy] = (float)(ix * (tamx - ix - 1) * iy * (tamy - iy - 1));
}

void
prtdat(tamx, tamy, u1)
long unsigned tamx, tamy;
float *u1;
{
   long unsigned i, j;
   // FILE *fp;

	// fp = fopen(fnam, "w");
	for (j = 0; j < tamx; ++j) {
		for (i = 0; i < tamy; ++i) {
			printf("%g", *(u1 + i*tamx + j) );
			// printf("i %lu j%lu\t", i, j);
			if (i != tamx-1){
				printf(" ");
			}else{
				printf("\n");
			}
		}
	}

	// fclose(fp);
	// printf(" %s\n",fnam);
}
