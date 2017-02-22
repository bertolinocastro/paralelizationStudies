#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>
#include <mpi.h>
#include <math.h>

#define C0 0.5
#define C1 0.75
#define C2 0.55
#define C3 0.25
#define C4 0.05
#define INICIO 1
#define FIM 2

// #define VERBOSE

int write_output(float *, int, int);
float get_time(int);

void printfcomma (long long unsigned n);

struct timeval inicio, final;

int main(int argc, char* argv[]){

	int niter,dim,num_th;
	int t,i,j,k,t0,t1; long long unsigned ini, fim, altura, resto;
	long double *heat_/*, tf_sec*/;

	MPI_Init( &argc, &argv );

		int rank, size;

		MPI_Comm_rank( MPI_COMM_WORLD, &rank );
		MPI_Comm_size( MPI_COMM_WORLD, &size );

		MPI_Request request1[size], request2[size];
		MPI_Status 	status1[size], 	status2[size];

		MPI_Group 	grupo_seleto;
		MPI_Comm 	comm_seleto;

		//executei /a.out 2 16 2
		if(argc < 4){
			printf("Error! informe o número de iterações e o tamanho da matriz. i.e: ./executavel 1000 500 500\n");
			MPI_Abort( MPI_COMM_WORLD, -1 );
		}
		niter 	= atoi(argv[1]);
		dim 	= atoi(argv[2]);
		num_th 	= atoi(argv[3]);

		omp_set_num_threads( num_th );

		altura = dim/size;

		/*
		 *
		 *	Verificar se a altura e menor que quatro, se sim, diminuir numero de processos.
		 *
		 *
		 */

		if( altura < 4 ){

			int h = size - 1;
			while( dim/h < 4 ) h--;

			printf("Numero de processos alto demais! Serao destruidos %d.\n", size - h);

			MPI_Comm_group(MPI_COMM_WORLD, &grupo_seleto);

			MPI_Group_range_excl(grupo_seleto, 1, &(int[3]){h, size-1,1}, &grupo_seleto);

			MPI_Comm_create(MPI_COMM_WORLD, grupo_seleto, &comm_seleto);

		}else{

			MPI_Comm_group(MPI_COMM_WORLD, &grupo_seleto);

			MPI_Comm_create(MPI_COMM_WORLD, grupo_seleto, &comm_seleto);


		}

		if( comm_seleto == MPI_COMM_NULL ){
			MPI_Finalize();
			exit(0);
		}

		MPI_Comm_rank( comm_seleto, &rank );
		MPI_Comm_size( comm_seleto, &size );

		altura = dim/size;
		resto = dim%size;


		/*
		*	Criacao do passo vertical: divisao da dimensao pelo size.
		*	Ex.: dim == 20 & size == 4 => alt == (5 p/proc.)+4
		*/

		if( size > 1 ){
			if(!rank){
				ini = 0;
				fim = 0 + altura + 4;
			}else if(rank == size-1){
				ini = 4;
				fim = 4 + altura + 0 + resto;
			}else{
				ini = 4;
				fim = 4 + altura + 4;
			}
		}else{
			ini = 0;
			fim = altura;
		}
		


		/*
		 *
		 *	Definicao dos processos vizinhos
		 *
		 */

		int left 	= (!rank)		 ? MPI_PROC_NULL : rank - 1,
			right 	= (rank==size-1) ? MPI_PROC_NULL : rank + 1;



		// altura_central = (dim/size);

		// altura_superior = altura_central-4;
		// altura_inferior = altura_central+4;

		//if(floor((float)dim/size)!=(float)dim/size){ inexato = 1; printf("VALOR DE ALTURA NAO DIVISIVEL. FALTA IMPLEMENTAR.\n"); MPI_Abort( MPI_COMM_WORLD, -1 ); }

		heat_ = (long double*) malloc(sizeof(long double)*2*fim*dim*dim);
		if(heat_ == NULL){ printf("Error! Malloc fail\n"); MPI_Abort( comm_seleto, -1 ); }

		long double (*heat)[fim][dim][dim] = (long double (*)[fim][dim][dim]) heat_;

		#pragma omp parallel for schedule(dynamic) private(i,j,k)
		for(i=0;i<fim;i++){
			for(j=0;j<dim;j++){
				for(k=0;k<dim;k++){
					heat[1][i][j][k] = heat[0][i][j][k] = 1.0f;
				}
			}
		}

		// if(!rank){

		// 	MPI_Isend(&heat[fim-4][0][0], 4*dim*dim, MPI_FLOAT, 1, 0, MPI_COMM_WORLD, request[rank]);
		// 	MPI_Recv(&heat[fim][0][0], 4*dim*dim, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, status[rank]);

		// }else if(rank==size-1){

		// 	MPI_Isend(&heat[ini][0][0], 4*dim*dim, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, request[rank]);
		// 	MPI_Recv(&heat[ini-4][0][0], 4*dim*dim, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, status[rank]);

		// }else{

		// 	//Enviando parcela propria ao antecessor
		// 	MPI_Isend(&heat[ini][0][0], 4*dim*dim, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, request[rank]);
		// 	//Enviando parcela propria ao sucessor
		// 	MPI_Isend(&heat[fim-4][0][0], 4*dim*dim, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD, request[rank]);

		// 	//Recebendo parcela do antecessor
		// 	MPI_Recv(&heat[ini-4][0][0], 4*dim*dim, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, request[rank]);
		// 	//Recebendo parcela do sucessor
		// 	MPI_Recv(&heat[fim][0][0], 4*dim*dim, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD, request[rank]);
		
		// }


/*
int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
             MPI_Comm comm)


int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm, MPI_Request *request)

int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
             MPI_Comm comm, MPI_Status *status)


*/

		printf("Sou o proc %d com ini %llu fim %llu altura %llu e resto %llu\n", rank, ini, fim, altura, resto),
		printf("Sou o proc %d e meu heat_ consome ", rank), printfcomma(sizeof(*heat)/8), printf("B de memoria\n");

		//printf("Oi, sou o processo %d e estou antes do loop principal!\n\n", rank);

		//get_time(INICIO);

		double start = MPI_Wtime();
		for(t=1; t<niter+1; t++){
			t0 = (t % 2);
			t1 = (t0 + 1)%(2);

			#pragma omp parallel for schedule(runtime) private(i,j,k)
			for(i=4;i<fim-4;i++){
			  for(j=4;j<dim-4;j++){
				for(k=4;k<dim-4;k++){
					heat[t0][i][j][k] = C0 * heat[t1][i][j][k]
								+ C1 * (heat[t1][i+1][j][k] + heat[t1][i-1][j][k] + heat[t1][i][j+1][k] + heat[t1][i][j-1][k] + heat[t1][i][j][k+1] + heat[t1][i][j][k-1])
								+ C2 * (heat[t1][i+2][j][k] + heat[t1][i-2][j][k] + heat[t1][i][j+2][k] + heat[t1][i][j-2][k] + heat[t1][i][j][k+2] + heat[t1][i][j][k-2])
								+ C3 * (heat[t1][i+3][j][k] + heat[t1][i-3][j][k] + heat[t1][i][j+3][k] + heat[t1][i][j-3][k] + heat[t1][i][j][k+3] + heat[t1][i][j][k-3])
								+ C4 * (heat[t1][i+4][j][k] + heat[t1][i-4][j][k] + heat[t1][i][j+4][k] + heat[t1][i][j-4][k] + heat[t1][i][j][k+4] + heat[t1][i][j][k-4]);
				}
			  }
			}

			MPI_Isend(&heat[t0][4][0][0], 4*dim*dim, MPI_LONG_DOUBLE, left, 0, comm_seleto, &request1[rank]);
			MPI_Isend(&heat[t0][ini+altura-4][0][0], 4*dim*dim, MPI_LONG_DOUBLE, right, 1, comm_seleto, &request2[rank]);


			MPI_Recv(&heat[t0][0][0][0], 4*dim*dim, MPI_LONG_DOUBLE, left, 1, comm_seleto, &status1[rank]);

			MPI_Recv(&heat[t0][ini+altura][0][0], 4*dim*dim, MPI_LONG_DOUBLE, right, 0, comm_seleto, &status2[rank]);

			/* RETIRAR O MPI_WAIT AUMENTA EM 30% O CUSTO TEMPORAL */
			MPI_Wait(&request1[rank], &status1[rank]);
			MPI_Wait(&request2[rank], &status2[rank]);

		}

		double end = MPI_Wtime();

		//tf_sec = get_time(FIM);


		// get_time(INICIO);
		// for(t=1; t<niter+1; t++){
		// 	t0 = (t % 2);
		// 	t1 = (t0 + 1)%(2);
			


		// 	#pragma omp parallel for schedule(runtime) private(i,j,k)
		// 	for(i=4;i<altura_inferior-4;i++){
		// 	  for(j=4;j<dim-4;j++){
		// 		for(k=4;k<dim-4;k++){
		// 			heat[t0][i][j][k] = C0 * heat[t1][i][j][k]
		// 						+ C1 * (heat[t1][i+1][j][k] + heat[t1][i-1][j][k] + heat[t1][i][j+1][k] + heat[t1][i][j-1][k] + heat[t1][i][j][k+1] + heat[t1][i][j][k-1])
		// 						+ C2 * (heat[t1][i+2][j][k] + heat[t1][i-2][j][k] + heat[t1][i][j+2][k] + heat[t1][i][j-2][k] + heat[t1][i][j][k+2] + heat[t1][i][j][k-2])
		// 						+ C3 * (heat[t1][i+3][j][k] + heat[t1][i-3][j][k] + heat[t1][i][j+3][k] + heat[t1][i][j-3][k] + heat[t1][i][j][k+3] + heat[t1][i][j][k-3])
		// 						+ C4 * (heat[t1][i+4][j][k] + heat[t1][i-4][j][k] + heat[t1][i][j+4][k] + heat[t1][i][j-4][k] + heat[t1][i][j][k+4] + heat[t1][i][j][k-4]);
		// 		}
		// 	  }
		// 	}
		// }


		// tf_sec = get_time(FIM);
		//printf("Time elapsed: %.10f seg\n",tf_sec);

		printf("Time elapsed: %.10f seg\n", end-start);
		MPI_Barrier( comm_seleto );

		if( size < 2 )
			// Printando o proprio
			for(unsigned int m=0;m<fim;m++){
				for(j=0;j<dim;j++){
					for(k=0;k<dim;k++){
						// printf("i %d j %d k %d heat %g\n",m,j,k,heat[t0][m][j][k]);
						printf("%LG\n",heat[t0][m][j][k]);

					}
				}
			}
		else if( !rank ){
			unsigned int m;
			// Printando o proprio
			for(m=0;m<fim-4;m++){
				for(j=0;j<dim;j++){
					for(k=0;k<dim;k++){
						// printf("i %d j %d k %d heat %g\n",m,j,k,heat[t0][m][j][k]);
						printf("%LG\n",heat[t0][m][j][k]);

					}
				}
			}

			heat_ = (long double*) realloc( heat_, sizeof(long double)*(altura+8+resto)*dim*dim );
			if( !heat_ ){ printf("Problema no realloc para impressao!\n"); MPI_Abort( MPI_COMM_WORLD, -1 ); }

			long double (*vetor)[altura+8+resto][dim][dim] = (long double (*)[altura+8+resto][dim][dim]) heat_;

			for( int l = 1; l < size; l++ ){
				unsigned int topo = ( l != size-1 ? altura+8 : altura + 4 + resto );
				unsigned int tipo = ( l != size-1 ? 4 : 0 );

				MPI_Recv( &vetor[0][0][0][0], dim*dim*topo, MPI_LONG_DOUBLE, l, 10+l, MPI_COMM_WORLD, &status1[0] );
				//printf("Recebido do rank %d topo vale: %u\n", l, topo);
				for(i=4;i<topo-tipo;i++,m++){
					for(j=0;j<dim;j++){
						for(k=0;k<dim;k++){
							// printf("i %d j %d k %d heat %g\n",m,j,k,vetor[t0][i][j][k]);
							printf("%LG\n",vetor[t0][i][j][k]);
						}
					}
				}

			}

		}else{
			MPI_Send( &heat[t0][0][0][0], fim*dim*dim, MPI_LONG_DOUBLE, 0, 10+rank, MPI_COMM_WORLD );
		}

		free(heat_);

	MPI_Finalize();

	return 0;
}
 
float get_time(int mode){

	float tf_sec = 0.0;

	if(mode == INICIO){
		gettimeofday(&inicio, NULL);
	}
	else{
		gettimeofday(&final, NULL);
		unsigned long long seg = 1000 * (final.tv_sec - inicio.tv_sec) + (final.tv_usec - inicio.tv_usec) / 1000;
		tf_sec = ((float)seg)*1e-3;
	}
	return tf_sec;
}

int write_output(float *heat_, int dim, int timestep){

	// FILE *fff;
	int i,j,k;
	float (*heat)[dim][dim][dim] = (float (*)[dim][dim][dim]) heat_;

	// fff = fopen("output.txt","w");
	// if (fff==NULL){
	// 	printf("Error! Não foi possível abrir o arquivo para saída.\n");
	// 	return 1;
	// }

	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++){
			for(k=0;k<dim;k++){
				printf("i %d j %d k %d heat %g\n",i,j,k,heat[timestep][i][j][k]);
			}
		}
	}
	// fclose(fff);
	return 0;
}

void printfcomma2 (int n) {
    if (n < 1000) {
        printf ("%d", n);
        return;
    }
    printfcomma2 (n/1000);
    printf (",%03d", n%1000);
}

void printfcomma (long long unsigned n) {
    if (n < 0) {
        printf ("-");
        n = -n;
    }
    printfcomma2 (n);
}
