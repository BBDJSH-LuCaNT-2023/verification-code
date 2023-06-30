all: phase1 phase2

phase1: phase1.c
	gcc -Iarcht -O -o phase1 phase1.c -lmpfi -lmpfr -lgmp -lm

phase2: phase2.c
	mpicc -Iarcht -O -o phase2 phase2_mpi.c -lmpfi -lmpfr -lgmp -lm

clean:
	rm -rf phase1 phase2
