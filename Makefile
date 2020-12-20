all:
	mpicxx -std=c++11 -O3 -fopenmp Parallel_Lab_2_MPI.cpp

clean:
	rm a.out
