#include <stdio.h>
#include <omp.h>
#include <stdlib.h>

int main(int argc, char const *argv[])
{
	int nth = strtol( argv[1], NULL, 10 );
	int x = 0;
	if(argc < 2) exit(1);
	omp_set_num_threads( nth );
	
	printf("Sou o x antes da regiao parallela %d\n", x);
	#pragma omp parallel
	{
		int id = omp_get_thread_num();
		int tot = omp_get_num_threads();

		#pragma omp sections reduction(+:x)
		{
			#pragma omp section
				printf("[1], id %d e do tot %d x vale: %d\n", id, tot, x=10);
			#pragma omp section
				printf("[2], id %d e do tot %d x vale: %d\n", id, tot, x=3);
			#pragma omp section
				printf("[3], id %d e do tot %d x vale: %d\n", id, tot, x=2);
			#pragma omp section
				printf("[4], id %d e do tot %d x vale: %d\n", id, tot, x=5);
		}

		printf("\tSou o th %d e meu x vale: %d\n", id, x);

	}
	printf("Sou o x apos a regiao parallela %d\n", x);


	return 0;
}