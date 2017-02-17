#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>
#include <mpi.h>

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

struct timeval inicio, final;

int main(int argc, char* argv[]){

	int niter,dim,num_th, altura;
	int t,i,j,k,t0,t1;
	float *heat_, tf_sec;

	MPI_Init( &argc, &argv );

		int rank, size;

		MPI_Comm_rank( MPI_COMM_WORLD, &rank );
		MPI_Comm_size( MPI_COMM_WORLD, &size );

		MPI_Request request[size];
		MPI_Status status[size];

		//executei /a.out 2 16 2
		if(argc < 4){
			printf("Error! informe o número de iterações e o tamanho da matriz. i.e: ./executavel 1000 500 500\n");
			exit(1);
		}
		niter 	= atoi(argv[1]);
		dim 	= atoi(argv[2]);
		num_th 	= atoi(argv[3]);

		omp_set_num_threads( num_th );

		/*
		*	Criacao do passo vertical: divisao da dimensao pelo size.
		*	Ex.: dim == 20 & size == 4 => alt == (5 p/proc.)+4
		*/

		if(!rank){
			inicio = 0;
			fim = altura + 4;
		}else if(rank == size-1){
			inicio = 4; // Equivalente a 0
			fim = inicio + altura;
		}else{
			inicio = 4; // equivalente a 0;
			fim = inicio + altura;
		}

		// altura_central = (dim/size);

		// altura_superior = altura_central-4;
		// altura_inferior = altura_central+4;

		if(floor((float)dim/size)!=(float)dim/size)	printf("VALOR DE ALTURA NAO DIVISIVEL. FALTA IMPLEMENTAR.\n");

		heat_ = (float*) malloc(sizeof(float)*2*(fim+4)*dim*dim);
		if(heat_ == NULL){ printf("Error! Malloc fail\n"); exit(1); }



		float (*heat)[altura_inferior][dim][dim] = (float (*)[altura_inferior][dim][dim]) heat_;

		#pragma omp parallel for schedule(dynamic) private(i,j,k)
		for(i=0;i<altura_inferior;i++){
			for(j=0;j<dim;j++){
				for(k=0;k<dim;k++){
					heat[1][i][j][k] = heat[0][i][j][k] = 1;
				}
			}
		}

		// if(!rank){

		// 	MPI_Isend(&heat[fim-4][0][0], 4*dim*dim, MPI_FLOAT, 1, 0, MPI_COMM_WORLD, request[rank]);
		// 	MPI_Recv(&heat[fim][0][0], 4*dim*dim, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, status[rank]);

		// }else if(rank==size-1){

		// 	MPI_Isend(&heat[inicio][0][0], 4*dim*dim, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, request[rank]);
		// 	MPI_Recv(&heat[inicio-4][0][0], 4*dim*dim, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, status[rank]);

		// }else{

		// 	//Enviando parcela propria ao antecessor
		// 	MPI_Isend(&heat[inicio][0][0], 4*dim*dim, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, request[rank]);
		// 	//Enviando parcela propria ao sucessor
		// 	MPI_Isend(&heat[fim-4][0][0], 4*dim*dim, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD, request[rank]);

		// 	//Recebendo parcela do antecessor
		// 	MPI_Recv(&heat[inicio-4][0][0], 4*dim*dim, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, request[rank]);
		// 	//Recebendo parcela do sucessor
		// 	MPI_Recv(&heat[fim][0][0], 4*dim*dim, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD, request[rank]);
		
		// }

		if( rank )
		//Enviando parcela propria ao antecessor
		MPI_Isend(&heat[inicio][0][0], 4*dim*dim, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, request[rank]);

		if( rank != size-1 )
		//Enviando parcela propria ao sucessor
		MPI_Isend(&heat[fim-4][0][0], 4*dim*dim, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD, request[rank]);

		if( rank )
		//Recebendo parcela do antecessor
		MPI_Recv(&heat[inicio-4][0][0], 4*dim*dim, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, request[rank]);
	
		if( rank != size-1 )
		//Recebendo parcela do sucessor
		MPI_Recv(&heat[fim][0][0], 4*dim*dim, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD, request[rank]);





/*
int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
             MPI_Comm comm)


int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm, MPI_Request *request)

int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
             MPI_Comm comm, MPI_Status *status)


*/

		get_time(INICIO);
		for(t=1; t<niter+1; t++){
			t0 = (t % 2);
			t1 = (t0 + 1)%(2);
			


			#pragma omp parallel for schedule(runtime) private(i,j,k)
			for(i=4;i<altura_inferior-4;i++){
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
		}
		tf_sec = get_time(FIM);





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

		// /*
		// *	Tentativa de segmentar o cubo em cubos menores com base no numero de threads
		// */




		// tf_sec = get_time(FIM);
		printf("Time elapsed: %.10f seg\n",tf_sec);
		//write_output( heat_, dim, t0 );

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
