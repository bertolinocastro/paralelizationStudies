#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define C0 0.5
#define C1 0.75
#define C2 0.55
#define C3 0.25
#define C4 0.05
#define INICIO 1
#define FIM 2

int write_output(long double *, int, int);
float get_time(int);

struct timeval inicio, final;

int main(int argc, char* argv[]){

	int niter,dim,save;
	int t,i,j,k,t0,t1;
	long double *heat_; float tf_sec;

	//executei /a.out 2 16 2
	if(argc < 4){
	    printf("Error! informe o número de iterações e o tamanho da matriz. i.e: ./executavel 1000 500 500\n");
	exit(1);
	}
	niter = atoi(argv[1]);
	dim = atoi(argv[2]);
	save = atoi(argv[3]);

	#ifdef VERBOSE
	printf("Número de iterações %d\ndimensão %d salvo a cada %d timesteps\n",niter,dim,save);
	#endif

	heat_ = (long double*) malloc(sizeof(long double)*2*dim*dim*dim);

	if(heat_ == NULL){
	    printf("Error! Malloc fail\n");
	    exit(1);
	}

	long double (*heat)[dim][dim][dim] = (long double (*)[dim][dim][dim]) heat_;

	#ifdef VERBOSE
	printf("Inicializando matriz!\n");
	#endif

	for(i=0;i<dim;i++){
	for(j=0;j<dim;j++){
	 for(k=0;k<dim;k++){
	heat[1][i][j][k] = heat[0][i][j][k] = 1;
	 }
	}
	}

	#ifdef VERBOSE
	printf("Iniciando computação do stencil!\n");
	#endif

	get_time(INICIO);
	for(t=1; t<niter+1; t++){
	    t0 = (t % 2);
	    t1 = (t0 + 1)%(2);
	    for(i=4;i<dim-4;i++){
	      for(j=4;j<dim-4;j++){
	        for(k=4;k<dim-4;k++){
	            heat[t0][i][j][k] = C0 * heat[t1][i][j][k] + C1 * (heat[t1][i+1][j][k] + heat[t1][i-1][j][k] + heat[t1][i][j+1][k] + heat[t1][i][j-1][k] + heat[t1][i][j][k+1] + heat[t1][i][j][k-1])
							   + C2 * (heat[t1][i+2][j][k] + heat[t1][i-2][j][k] + heat[t1][i][j+2][k] + heat[t1][i][j-2][k] + heat[t1][i][j][k+2] + heat[t1][i][j][k-2])
							   + C3 * (heat[t1][i+3][j][k] + heat[t1][i-3][j][k] + heat[t1][i][j+3][k] + heat[t1][i][j-3][k] + heat[t1][i][j][k+3] + heat[t1][i][j][k-3])
							   + C4 * (heat[t1][i+4][j][k] + heat[t1][i-4][j][k] + heat[t1][i][j+4][k] + heat[t1][i][j-4][k] + heat[t1][i][j][k+4] + heat[t1][i][j][k-4]);
	        }
	      }
	    }
	}
	tf_sec = get_time(FIM);
	printf("Time elapsed: %.10f seg\n",tf_sec);
	write_output(heat_, dim, t0);

	free(heat_);

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

int write_output(long double *heat_, int dim, int timestep){
 
 // FILE *fff;
 int i,j,k;
 long double (*heat)[dim][dim][dim] = (long double (*)[dim][dim][dim]) heat_;
 
 // fff = fopen("output.txt","w");
 // if (fff==NULL){
 //   printf("Error! Não foi possível abrir o arquivo para saída.\n");
 //   return 1;
 // }

 for(i=0;i<dim;i++){
   for(j=0;j<dim;j++){
     for(k=0;k<dim;k++){
		// printf("i %d j %d k %d heat %g\n",i,j,k,heat[timestep][i][j][k]);
		printf("%LG\n",heat[timestep][i][j][k]);
	}
   }
 }
 // fclose(fff);
 return 0;
}
