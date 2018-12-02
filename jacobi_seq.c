/**
* Author: Ruben Gran Tejero
* 
* Description: Sequential version of a simulation of thermal transmission in a 2D space.
*
* If you use this code for academic work, please also reference:
*      Ruben Gran Tejero, rgran@unizar.es.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>

//Reference implementation of Jacobi for mp2 (sequential)
//Constants are being used instead of arguments
#define BC_HOT  1.0
#define BC_COLD 0.0
#define INITIAL_GRID 0.5
#define N_DIM 16 
#define MAX_ITERATIONS 1000
#define TOL 1.0e-4

struct timeval tv;
double get_clock() {
   struct timeval tv; int ok;
   ok = gettimeofday(&tv, (void *) 0);
   if (ok<0) { printf("gettimeofday error");  }
   return (tv.tv_sec * 1.0 + tv.tv_usec * 1.0E-6);
}


double **create_matrix(int n) {
	int i;
	double **a;

	a = (double**) malloc(sizeof(double*)*n);
	for (i=0;i<n;i++) {
		a[i] = (double*) malloc(sizeof(double)*n);
	}

	return a;
}

void init_matrix(double **a, int n) {

	int i, j;
	
	for(i=0; i<n; i++) {
		for(j=0; j<n; j++)
			a[i][j] = INITIAL_GRID;
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
			printf("%6.4lf ",a[i][j]);
		}
		printf("\n");
	}
}

void free_matrix(double **a, int n) {
	int i;
	for (i=0;i<n;i++) {
		free(a[i]);
	}
	free(a);
}

int main(int argc, char* argv[]) {
	int i,j,iteration;
	int n = N_DIM;
	double **a, **b, maxdiff;
	double tstart, tend, ttotal;

	//add 2 to each dimension to use sentinal values
	a = create_matrix(n+2);
	b = create_matrix(n+2);

	init_matrix(a,n+2);

	//Initialize the hot boundaries
	for(i=0;i<n+2;i++) {
		a[i][0] = BC_HOT;
	    a[i][n+1] = BC_HOT;
	    a[0][i] = BC_HOT;
	}

	// Initialize the cold boundary
	for(j=0;j<n+2;j++) {
		a[n+1][j] = BC_COLD;
	}

	// Copy a to b
	for(i=0; i<n+2; i++) {
		for(j=0; j<n+2; j++) {
			b[i][j] = a[i][j];
		}
	}

	// Main simulation routine
	iteration=0;
	maxdiff=1.0;
	printf("Running simulation with tolerance=%lf and max iterations=%d\n",
		TOL, MAX_ITERATIONS);
	tstart = get_clock();
	int iter=0;
	while(maxdiff > TOL && iteration<MAX_ITERATIONS) {

		//printf("Iteration=%d\n",iter);
		iter++;

		// Compute new grid values
		maxdiff = 0.0;
		for(i=1;i<n+1;i++) {
			for(j=1;j<n+1;j++) {
				b[i][j] = 0.2*(a[i][j]+a[i-1][j]+a[i+1][j]+a[i][j-1]+a[i][j+1]);
		        if (fabs(b[i][j]-a[i][j]) > maxdiff)
		          maxdiff = fabs(b[i][j]-a[i][j]);
			}
		}

		// Copy b to a
		swap_matrix(&a,&b);	

		iteration++;
	}
	tend = get_clock();
	ttotal = tend-tstart;

	// Output final grid
	printf("Final grid:\n");
	print_grid(a,0,n+2);

	// Results
	printf("Results:\n");
	printf("Iterations=%d\n",iteration);
	printf("Tolerance=%12.10lf\n",maxdiff);
	printf("Running time=%12.10lf\n",ttotal);

	free_matrix(a,n+2);
	free_matrix(b,n+2);
	return 0;
}
