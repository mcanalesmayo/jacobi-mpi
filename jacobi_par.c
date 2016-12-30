/**
* Author: Marcos Canales Mayo
* 
* Description: MPI (parallel) version of jacobi_seq.c
*
* If you use this code for academic work, please also reference:
*      Ruben Gran Tejero, rgran@unizar.es.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <mpi.h>
#include <string.h>
#include <float.h>

// Problem constants
#define BC_HOT  1.0
#define BC_COLD 0.0
#define INITIAL_GRID 0.5
// Jacobi constants
#define MAX_ITERATIONS 1000
#define TOL 1.0e-4

struct timeval tv;
double get_clock() {
   struct timeval tv; int ok;
   ok = gettimeofday(&tv, (void *) 0);
   if (ok<0) { printf("gettimeofday error");  }
   return (tv.tv_sec * 1.0 + tv.tv_usec * 1.0E-6);
}

double **create_matrix(int subprob_size) {
	int i;
	double **a;
	double *rows;

	// in this approach matrix is stored as array of pointers to array of doubles
	a = (double**) malloc(sizeof(double*)*subprob_size);
	// ensure the matrix is contiguous in memory so it can be accessed as a static matrix: a[i][j]
	// allocating the matrix in a loop, row by row, may not work as reserved memory is not ensured to be contiguous
	rows = (double*) malloc(sizeof(double)*subprob_size*subprob_size);

	for (i=0;i<subprob_size;i++) {
		a[i] = &rows[i*subprob_size];
	}

	return a;
}

void init_matrix(double **a, double *rfrbuff, double *rfcbuff, double *rlrbuff, double *rlcbuff, int n_subprobs, int subprob_size, int column_num, int row_num) {
	int i, j;

	// Initialize matrix
	// First time all values are INITIAL_GRID
	for(i=0; i<subprob_size; i++) {
		rfrbuff[i] = INITIAL_GRID;
		rfcbuff[i] = INITIAL_GRID;
		rlrbuff[i] = INITIAL_GRID;
		rlcbuff[i] = INITIAL_GRID;
		for(j=0; j<subprob_size; j++)
			a[i][j] = INITIAL_GRID;
	}

	// Switch requires compile time knowledge of cases
	// Alternative is if
	// Outline column values
	// I'm in the first column
	if (column_num == 0){
		for(i=0; i<subprob_size; i++){
			rfcbuff[i] = BC_HOT;
		}
	}
	// I'm in the last column
	else if(column_num == ((int) sqrt(n_subprobs))-1){
		for(i=0; i<subprob_size; i++){
			rlcbuff[i] = BC_HOT;
		}
	}

	// Outline row values
	// I'm in the first row
	if (row_num == 0){
		for(j=0; j<subprob_size; j++){
			rfrbuff[j] = BC_HOT;
		}
	}
	// I'm in the last row
	else if(row_num == ((int) sqrt(n_subprobs))-1){
		for(j=0; j<subprob_size; j++){
			rlrbuff[j] = BC_COLD;
		}
	}
}

void swap_matrix(double ***a, double ***b) {
	double **temp;

	temp = *a;
	*a = *b;
	*b = temp;
}

void print_grid(double **a, int nstart, int nend) {
	int i, j;

	for(i=nstart; i<nend; i++) {
		for(j=nstart; j<nend; j++) {
			printf("%6.4lf ", a[i][j]);
		}
		printf("\n");
	}
}

void free_matrix(double **a) {
	int i;

	// ** reverse order would raise segm. fault
	// ** as pointers to rows would be freed before rows themselves are freed
	free(a[0]);
	free(a);
}

void wait_req(MPI_Request *req){
	int req_flag = 0;
	while(!req_flag) MPI_Test(req, &req_flag, MPI_STATUS_IGNORE);
}

int main(int argc, char* argv[]) {
	// instead of using always the same tag, in order to provide more information about messages
	// if needed 'iteration' is also used as MPI tag
	int i, j , i_aux = 0, j_aux = 0, generic_tag = 0, iteration;
	int n_dim, n_subprobs, subprob_size;
	int column_num, row_num;
	double **a, **b, maxdiff, maxdiff_aux;
	MPI_Datatype double_strided_vect;
	// root process
	double **res;
	int root_rank = 0, res_offset;

	double tstart, tend, ttotal;

	// MPI vars
	int my_rank;
	// [s]end/[r]eceive [f]irst/[l]ast [r]ow/[c]olumn [buff]ers
	double *sfrbuff, *sfcbuff, *slrbuff, *slcbuff;
	double *rfrbuff, *rfcbuff, *rlrbuff, *rlcbuff;
	// [s]end/[r]eceive [f]irst/[l]ast [r]ow/[c]olumn [req]uests
	MPI_Request *sfrreq, *sfcreq, *slrreq, *slcreq;
	MPI_Request *rfrreq, *rfcreq, *rlrreq, *rlcreq;

	if (argc != 3) return -1;

	// Init MPI lib
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	// My subproblem params
	n_subprobs = atoi(argv[1]);
	n_dim = atoi(argv[2]);
	
	subprob_size = (int) sqrt((n_dim*n_dim)/n_subprobs);

	// Where am I
	column_num = my_rank%((int) sqrt(n_subprobs));
	row_num = (int) (my_rank/(int) sqrt(n_subprobs));

	// Could do a single big malloc to avoid overhead of multiple syscalls
	sfrbuff = malloc(subprob_size*sizeof(double)); sfcbuff = malloc(subprob_size*sizeof(double)); slrbuff = malloc(subprob_size*sizeof(double)); slcbuff = malloc(subprob_size*sizeof(double));
	rfrbuff = malloc(subprob_size*sizeof(double)); rfcbuff = malloc(subprob_size*sizeof(double)); rlrbuff = malloc(subprob_size*sizeof(double)); rlcbuff = malloc(subprob_size*sizeof(double));

	sfrreq = malloc(sizeof(MPI_Request)); sfcreq = malloc(sizeof(MPI_Request)); slrreq = malloc(sizeof(MPI_Request)); slcreq = malloc(sizeof(MPI_Request));
	rfrreq = malloc(sizeof(MPI_Request)); rfcreq = malloc(sizeof(MPI_Request)); rlrreq = malloc(sizeof(MPI_Request)); rlcreq = malloc(sizeof(MPI_Request));

	// Alloc matrices
	a = create_matrix(subprob_size);
	b = create_matrix(subprob_size);
	if (my_rank == root_rank) res = create_matrix(n_dim);

	// Create strided vector datatype, used when gathering all subproblems
	MPI_Type_vector(subprob_size, subprob_size, n_dim, MPI_DOUBLE, &double_strided_vect);
	MPI_Type_commit(&double_strided_vect);

	// Main simulation routine
	iteration=0;
	printf("[%d] Running simulation with tolerance=%lf and max iterations=%d\n",
		my_rank, TOL, MAX_ITERATIONS);
	tstart = MPI_Wtime();

	init_matrix(a, rfrbuff, rfcbuff, rlrbuff, rlcbuff, n_subprobs, subprob_size, column_num, row_num);

	maxdiff = DBL_MAX;
	while(maxdiff > TOL && iteration<MAX_ITERATIONS) {
		maxdiff = 0.0;

		// Send my outer columns values
		// I'm not in the last column
		if (column_num != ((int) sqrt(n_subprobs))-1){
			// send the last column of my subproblem matrix
			for(i=0; i<subprob_size; i++) slcbuff[i] = a[i][subprob_size-1];
			MPI_Isend(slcbuff, subprob_size, MPI_DOUBLE, my_rank+1, iteration, MPI_COMM_WORLD, slcreq);

			// post receive last column + 1
			MPI_Irecv(rlcbuff, subprob_size, MPI_DOUBLE, my_rank+1, iteration, MPI_COMM_WORLD, rlcreq);
		}
		// I'm not in the first column
		if (column_num != 0){
			// send the first column of my subproblem matrix
			for(i=0; i<subprob_size; i++) sfcbuff[i] = a[i][0];
			MPI_Isend(sfcbuff, subprob_size, MPI_DOUBLE, my_rank-1, iteration, MPI_COMM_WORLD, sfcreq);

			// post receive first column - 1
			MPI_Irecv(rfcbuff, subprob_size, MPI_DOUBLE, my_rank-1, iteration, MPI_COMM_WORLD, rfcreq);
		}

		// Send my outer rows values
		// I'm not in the last row
		if (row_num != ((int) sqrt(n_subprobs))-1){
			// send the last row of my subproblem matrix
			memcpy(slrbuff, a[subprob_size-1], subprob_size*sizeof(double));
			MPI_Isend(slrbuff, subprob_size, MPI_DOUBLE, (int) (my_rank+sqrt(n_subprobs)), iteration, MPI_COMM_WORLD, slrreq);

			// post receive last row + 1
			MPI_Irecv(rlrbuff, subprob_size, MPI_DOUBLE, (int) (my_rank+sqrt(n_subprobs)), iteration, MPI_COMM_WORLD, rlrreq);
		}
		// I'm not in the first row
		if (row_num != 0){
			// send the first row of my subproblem matrix
			memcpy(sfrbuff, a[0], subprob_size*sizeof(double));
			MPI_Isend(sfrbuff, subprob_size, MPI_DOUBLE, (int) (my_rank-sqrt(n_subprobs)), iteration, MPI_COMM_WORLD, sfrreq);

			// post receive first row - 1
			MPI_Irecv(rfrbuff, subprob_size, MPI_DOUBLE, (int) (my_rank-sqrt(n_subprobs)), iteration, MPI_COMM_WORLD, rfrreq);
		}

		// Compute new inner grid values
		// EL = 0.2*(   EL  +   UP     +  DOWN +   LEFT   + RIGHT );
		// Inner rows i=[1...subprob_size-2]
		for(i=1;i<subprob_size-1;i++) {
			// j=[1...subprob_size-2]
			for(j=1;j<subprob_size-1;j++) {
				b[i][j] = 0.2*(a[i][j]+a[i-1][j]+a[i+1][j]+a[i][j-1]+a[i][j+1]);
				if (fabs(b[i][j]-a[i][j]) > maxdiff) maxdiff = fabs(b[i][j]-a[i][j]);
			}
		}

		// Recv outline columns values
		// I'm not in the last column
		if (column_num != ((int) sqrt(n_subprobs))-1){
			// waiting to receive last column
			wait_req(rlcreq);
		}
		// I'm not in the first column
		if (column_num != 0){
			// waiting to receive first column
			wait_req(rfcreq);
		}

		// Recv outline rows values
		// I'm not in the last row row
		if (row_num != ((int) sqrt(n_subprobs))-1){
			// waiting to receive last row
			wait_req(rlrreq);
		}
		// I'm not in the first row
		if (row_num != 0){
			// waiting to receive first row
			wait_req(rfrreq);
		}

		// Compute new outer grid values
		// EL	= 0.2*(   EL  +   UP     +  DOWN +   LEFT   + RIGHT );
		// First row i=0
		i=0;
		// j=0
		j=0;
		b[i][j] = 0.2*(a[i][j]+rfrbuff[j]+a[i+1][j]+rfcbuff[i]+a[i][j+1]);
		if (fabs(b[i][j]-a[i][j]) > maxdiff) maxdiff = fabs(b[i][j]-a[i][j]);
		// j=[1...subprob_size-2]
		for(j=1;j<subprob_size-1;j++){
			b[i][j] = 0.2*(a[i][j]+rfrbuff[j]+a[i+1][j]+a[i][j-1]+a[i][j+1]);
			if (fabs(b[i][j]-a[i][j]) > maxdiff) maxdiff = fabs(b[i][j]-a[i][j]);
		}
		// j=subprob_size-1
		j=subprob_size-1;
		b[i][j] = 0.2*(a[i][j]+rfrbuff[j]+a[i+1][j]+a[i][j-1]+rlcbuff[i]);
		if (fabs(b[i][j]-a[i][j]) > maxdiff) maxdiff = fabs(b[i][j]-a[i][j]);

		// Inner rows i=[1...subprob_size-2]
		for(i=1;i<subprob_size-1;i++) {
			// j=0
			j=0;
			b[i][j] = 0.2*(a[i][j]+a[i-1][j]+a[i+1][j]+rfcbuff[i]+a[i][j+1]);
			if (fabs(b[i][j]-a[i][j]) > maxdiff) maxdiff = fabs(b[i][j]-a[i][j]);

			// j=subprob_size-1
			j=subprob_size-1;
			b[i][j] = 0.2*(a[i][j]+a[i-1][j]+a[i+1][j]+a[i][j-1]+rlcbuff[i]);
			if (fabs(b[i][j]-a[i][j]) > maxdiff) maxdiff = fabs(b[i][j]-a[i][j]);
		}

		// Last row i=subprob_size-1
		i=subprob_size-1;
		// j=0
		j=0;
		b[i][j] = 0.2*(a[i][j]+a[i-1][j]+rlrbuff[j]+rfcbuff[i]+a[i][j+1]);
		if (fabs(b[i][j]-a[i][j]) > maxdiff) maxdiff = fabs(b[i][j]-a[i][j]);
		// j=[1...subprob_size-2]
		for(j=1;j<subprob_size-1;j++){
			b[i][j] = 0.2*(a[i][j]+a[i-1][j]+rlrbuff[j]+a[i][j-1]+a[i][j+1]);
			if (fabs(b[i][j]-a[i][j]) > maxdiff) maxdiff = fabs(b[i][j]-a[i][j]);
		}
		// j=subprob_size-1
		j=subprob_size-1;
		b[i][j] = 0.2*(a[i][j]+a[i-1][j]+rlrbuff[j]+a[i][j-1]+rlcbuff[i]);
		if (fabs(b[i][j]-a[i][j]) > maxdiff) maxdiff = fabs(b[i][j]-a[i][j]);

		// All get the maximum diff
		MPI_Allreduce(&maxdiff, &maxdiff_aux, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		maxdiff = maxdiff_aux;

		// Copy b to a
		swap_matrix(&a,&b);

		iteration+=1;
	}

	// Gatherv doesn't fit here, as displs param isn't taken as bytes, instead it's multiplied by the size of the recvtype, which is double_strided_vect
	// Alternative is Send subprob_size*subprob_size contiguous doubles and Recv strided vector

	// Root process groups results
	if (my_rank == root_rank){
		// reuse sfrreq as req handler
		// send my values to myself -> IMPORTANT: non blocking mode
		MPI_Isend(a[0], subprob_size*subprob_size, MPI_DOUBLE, root_rank, generic_tag, MPI_COMM_WORLD, sfrreq);
		MPI_Recv(res[0], 1, double_strided_vect, root_rank, generic_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		for(i=1;i<n_subprobs;i++){
			// res_offset:	- add (row number * row size)
			// 				- add (in-row position * subproblem size)
			res_offset = ((int) (i/(int) sqrt(n_subprobs))) *subprob_size*subprob_size*(int) sqrt(n_subprobs);
			res_offset += i%((int) sqrt(n_subprobs)) * subprob_size;
			MPI_Recv(res[0] + res_offset, 1, double_strided_vect, i, generic_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}
	// Remaining processes
	else{
		// Send my subproblem results
		MPI_Send(a[0], subprob_size*subprob_size, MPI_DOUBLE, root_rank, generic_tag, MPI_COMM_WORLD);
	}

	tend = MPI_Wtime();
	ttotal = tend-tstart;

	if (my_rank == root_rank){
		// Output final grid
		printf("Final grid:\n");
		print_grid(res, 0, n_dim);

		// Results
		printf("Results:\n");
		printf("Iterations=%d\n", iteration);
		printf("Tolerance=%12.10lf\n", maxdiff);
		printf("Problem dimmensions=%dx%d\n", n_dim, n_dim);
		printf("Number of subproblems=%d\n", n_subprobs);
		printf("Running time=%12.10lf\n", ttotal);

		// Free allocated mem
		free_matrix(res);
	}
	free_matrix(a);
	free_matrix(b);

	MPI_Type_free(&double_strided_vect);

	// Finalize MPI lib
	MPI_Finalize();

	return 0;
}
