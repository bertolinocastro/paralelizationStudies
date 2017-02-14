CC=mpicc
CFLAGS=-fopenmp -Wall -pedantic -ansi -g -std=gnu99
LIBS=-lm

file: $(f).c
	$(CC) $(CFLAGS) $^ $(LIBS) -o $(f).out

clean:
	rm ./a.out