#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "../library/mmio.h"
#include "omp.h"
#include "../library/helper.h"

/**
* @description      this funtion is used for getting the row r with maximum value on column k locally, and return r if succeeded otherwise -1
* @param    a   	matrix A
* @param    k    	target column k
* @param    N 		the size of matrix
*/
int max_col(double *a, int k, int N){
	double max=fabs(*(a + k*N + k));//a_kk
	int i=0;
	int r=-1;
 
	for ( i = k+1; i < N; i++)
	{
		if (fabs(*(a + i* N + k)) >max ) // |a_ik|
		{
			r=i;
			max=fabs(*(a + i*N + k));
		}
	}
	return r;
}

/**
* @description      this funtion is used for exchanging row r and k in matrix A and b
* @param    a   	matrix A
* @param    b    	matrix b
* @param    r     	row r
* @param    k    	row k
* @param    N 		the size of matrix
*/
void exchange_row(double *a, double *b, int r, int k, int N){
	double temp;
	int i;
	for (i = 0; i < N; ++i)
	{
		temp=*(a+k*N+i);
		*(a+k*N+i)=*(a+r*N+i);
		*(a+r*N+i)=temp;
	}
	temp=*(b+k);
	*(b+k)=*(b+r);
	*(b+r)=temp;
}

/**
* @description      this funtion executes Gaussian elimination with shared memory using OpenMp
* @param    a   	matrix A
* @param    b    	matrix b
* @param    N 		the size of matrix
*/
double *gauss_openmp(double *a, double *b, int N){
	double *x;
	int k,r;
	int id , numthreads;
	x = (double *)malloc(N*sizeof(double));

	/* Forward elimination*/
	for ( k = 0; k < N-1; ++k)
	{
		r = max_col(a,k,N);

		if (r!=k && r!=-1)
		{
			exchange_row(a, b, r, k, N);
		}

		/* using OpenMp to execute elimination */
		#pragma omp parallel 
		{
			int id = omp_get_thread_num();
            int numthreads = omp_get_num_threads();
            for (int i = k+1+id; i < N; i+=numthreads)
            {
                // l[i] = a[i][k] / a[k][k];
                for (int j = k+1; j < N; j++)
                {
                	// a[i][j] = a[i][j] - l[i] * a[k][j];
                	*(a+i*N+j) = *(a+i*N+j) - ( *(a+i*N+k) / *(a+k*N+k) )*( *(a+k*N+j));
                }
                // b[i] = b[i] -l[i] * b[k];
                *(b+i) = *(b+i) - ( *(a+i*N+k) / *(a+k*N+k) )*(*(b+k)); 
            }
		}
	}
	double sum;
	int j;
	
	/* Backward substitution */
	for ( k = N-1; k>=0; k--)
	{
		sum=0;
		for (j = k+1; j < N; ++j)
		{
			sum += *(a+k*N+j) * *(x+j);
		}
		*(x+k) = 1/(*(a+k*N+k))*( *(b+k)-sum );
	}
	return x;
}

int main(int argc, char *argv[])
{
	int N = 8;//default
    double t1, t2;
    double start_time, run_time;

    int p1=1, p2=1;
    char *fname_in;
    char *fname_out;
    /* get command line options, set input and output file names and the value of p1 and p2 */
    handle_opt( argc, argv, &fname_in, &fname_out, &p1, &p2);

    double *A;
	/* read the matrix to A and the size of matrix size to N, then broadcast N to other processors*/
    A=input_matrix(fname_in, &N);

    double *b;
	b = (double *)malloc(N*sizeof(double));
	int i,j;
	/* all elements in matrix b are set to 1 */
    for ( i = 0; i < N; ++i)
	{
		*(b+i)=1;
	}

	start_time = omp_get_wtime();
	/* Executing Gaussian elimination with shared memory using OpenMp */
    double *x=gauss_openmp(A, b, N);

    run_time = omp_get_wtime() - start_time;

    /* write result x into output file and print the execution time */
	output_matrix(fname_out, x, N);
    printf("gauss_openmp() computation in %f seconds\n",run_time);
    return 0;
}
