#!/bin/bash

if [[ -z "$1" ]]; then
	echo
	echo "Script de diferenciacao de outputs"
	echo "Usabilidade: [<no. procs.> <path/file1>] [<no. procs.> <path/file2>]"
	echo
else
	echo "diffout:"
	diff <(mpiexec -np $1 $2) <(mpiexec -np $3 $4)
fi



