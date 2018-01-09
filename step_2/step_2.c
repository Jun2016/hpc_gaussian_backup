#include <stdlib.h>
#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include "mmio.h"
#include "helper.h"

#define MASTER 0 //set master rank to be rank 0
int N;

/**
* @description      this funtion is used for print A matrix
* @param    a    	matrix A
*/
void printA(double *a){
	int i,j;
	printf("A:\n");
	for ( i = 0; i < N; ++i)
	{
		for ( j = 0; j < N; ++j)
		{
			printf("%lf\t",*(a + i*N + j) );
		}
		printf("\n");
	}
	printf("\n");
}

/**
* @description      this funtion is used for getting the row r with maximum value on column k locally, and return r if succeeded otherwise -1
* @param    a   	matrix A
* @param    k    	target column k
* @param    p     	the number of processors
* @param    me    	the rank of current processor
* @param    loc_n   the number of rows belongs to current processor
*/
int max_col_loc(double *a, int k, int p,int me, int loc_n){
	double max=0;
	int i=0;
	int r=-1;

	while( me + i*p < k){
		i++;
	}
	int pos=me+i*p;
	
	if (pos>=N)
	{
		return -1;
	}
 
 	/* get the r which maximize |a_rk| */
	for ( i = pos; i < N; i+=p)
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
*/
void exchange_row(double *a, double *b, int r, int k){
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
* @description      this funtion is used for copying row k in matrix A and b into buf
* @param    a   	matrix A
* @param    b    	matrix b
* @param    k    	row k
* @param    buf     buffer
*/
void copy_row(double  *a, double *b, int k, double *buf){
	int i;
	for ( i = 0; i < N; ++i)
	{
		*(buf+i)=*(a + k*N +i);
	}
	*(buf+N)=*(b+k);
}

/**
* @description      this funtion is used for exchanging row r in matrix A and b with buffer buf
* @param    a   	matrix A
* @param    b    	matrix b
* @param    r    	row r
* @param    buf     buffer
*/
void copy_exchange_row(double *a, double *b, int r, double *buf){
	int i;
	double temp;

	for ( i = 0; i < N; ++i)
	{
		temp = *(a + r*N +i);
		*(a+r*N+i) = *(buf+i);
		*(buf+i) = temp;
	}
	temp = *(b+r);
	*(b+r) = *(buf+N);
	*(buf+N) = temp;
}

/**
* @description      this funtion is used for writting buf into row k in matrix A and b 
* @param    a   	matrix A
* @param    b    	matrix b
* @param    buf     buffer
* @param    k    	row k
*/
void copy_back_row(double *a, double *b, double *buf, int k){
	int i;
	for ( i = 0; i < N; ++i)
	{
		*(a + k*N +i) = *(buf+i);
	}
	*(b+k)=*(buf+N);
}

/**
* @description      this funtion executes Gaussian elimination with row-cyclic data distribution 
* @param    a   	matrix A
* @param    b    	matrix b
* @param    p    	the number of processors
* @param    me    	the rank of current processor
*/
double *gauss_cyclic(double *a, double *b, int p, int me){
	double *x, l[N], *buf;
	int i,j,k,r, tag=42;
	MPI_Status status;
	int loc_n=N/p; 

	struct{
		double val; 
		int node;
	}z,y;

	x = (double *)malloc(N*sizeof(double));
	buf = (double *) malloc((N+1)*sizeof(double));


	/* Forward elimination*/
	for ( k = 0; k < N-1; ++k)
	{
		/* get the local pivot row */
		r = max_col_loc(a, k, p, me, loc_n);
		/* put the maximum |a_rk| with rank of current processor into z */
		z.node = me;
		if (r != -1)
		{
			z.val = fabs(*(a + r*N + k));//|a_rk|
		}
		else{
			z.val= 0.0;
		}
		/* get the maximum |a_rk| with its processor's rank in y */
		MPI_Allreduce(&z, &y, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

		if(k%p==y.node){
			/* pivot row and row k are on the same processor */
			/* exchange row r and k and put the pivot row into buf */
			if (k % p == me)
			{
				if (*(a + k*N + k) != y.val) //a_kk
				{		
					exchange_row(a,b,r,k);
				}
				copy_row(a,b,k,buf);
			}
		}
		else{
			/* pivot row and row k are owned by different processors */
			if (k%p == me)
			{
				/* current processor has the row k */
				/* copy row k into buf and send to the processor has row r */
				copy_row(a,b,k,buf);
				MPI_Send(buf+k, N-k+1, MPI_DOUBLE, y.node, tag, MPI_COMM_WORLD);
			}
			else if(y.node == me){
				/* current processor has the row k */
				/* receive row k and exchange row r with buf */
				MPI_Recv(buf+k, N-k+1, MPI_DOUBLE, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
				copy_exchange_row(a,b,r,buf);
			}
		}
		/* broadcast the pivot row */
		MPI_Bcast(buf+k, N-k+1,MPI_DOUBLE, y.node, MPI_COMM_WORLD);

		/* if processor has row k but doesn't have row r, writting received pivot row into row k */
		if ((k%p != y.node) && (k%p==me))
		{
			copy_back_row(a,b,buf,k);
		}

		/* Eliminate A and b */
		i=k+1;
		while(i%p!=me)
			i++;
		for (; i < N; i+=p)
		{
			/* l_i = a_ik/buf_k */
			*(l+i)=(*(a+i*N+k))/(*(buf+k));
			for (j=k; j < N; ++j)
			{
				/* a_ij = a_ij - (l_i * buf_j) */
				*(a+i*N+j) = *(a+i*N+j) - (*(l+i)) * (*(buf+j));
			}
			/* b_i = b_i - (l_i * buf_N) */
			*(b+i) = *(b+i) - (*(l+i)) * (*(buf+N));
		}
	}
	
	double sum;
	for ( i = 0; i < N; ++i)
	{
			*(x+i)=0.;
	}

	/* Backward substitution */
	for (k=N-1; k>=0; k--)
	{
		/* If current processor has row k, calculate x_k*/
		if(k%p==me){
			sum=0.;
			for (j=k+1; j<N; ++j)
			{
				sum += *(a+k*N+j) * *(x+j);
			}
			*(x+k) = 1 / *(a+ k*N + k) * (*(b+k) - sum);
		}
		/* broadcast x_k among all processors*/
		MPI_Bcast(x+k, 1, MPI_DOUBLE, k%p, MPI_COMM_WORLD);
	}
	return x;
}


int main(int argc, char *argv[])
{
    int i, j, cur_msg_len, p, myrank;
    double t1, t2; 

    int p1=1, p2=1;
    char *fname_in;
    char *fname_out;

    /* get command line options, set input and output file names and the value of p1 and p2 */
    handle_opt( argc, argv, &fname_in, &fname_out, &p1, &p2);

    /* Initialize the MPI and get the number of processors into p */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    int n=0;

    /*if the distribution of processors is not specified, let it be row-cyclic districution */
    if (p1==-1)
    {
    	p1=p;
    	p2=1;
    }

    /* In step 2, p1 should be p and p2 should be 1 cause of row-cyclic districution */
    if (p2 != 1)
    {
		if(myrank==MASTER)printf("p2: %d will be set to 1\n", p2);
		p2 = 1;
    }
    if (p1 != p)
    {
    	if(myrank==MASTER)printf("ERROR: p1 not equals to p\n");
    	return 1;
    }

    /* get and put the rank of current processor to myrank */
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    double *A;

    /* only master rank processor read the matrix to A and the size of matrix size to N, then broadcast N to other processors*/
    if (myrank==MASTER)
    {
    	A=input_matrix(fname_in, &n);
    	N = n;
    }
    MPI_Bcast(&N, 1, MPI_INT, MASTER, MPI_COMM_WORLD); 
    MPI_Barrier(MPI_COMM_WORLD); 

    if (myrank!= MASTER)
	{
		A = (double *)malloc(N*N*sizeof(double));
	}

	/* all elements in matrix b are set to 1 */
	double *b;
	b = (double *)malloc(N*sizeof(double));
    for ( i = 0; i < N; ++i)
	{
		*(b+i)=1;
	}
    
    /*calculate the number of rows each processor holds*/
    int loc_n=N/p;

    /* Master rank processor broadcast matrix A to all the other processors */
	MPI_Bcast(A, N * N , MPI_DOUBLE, MASTER, MPI_COMM_WORLD);   
    MPI_Barrier(MPI_COMM_WORLD);

    t1 = MPI_Wtime();
    /* Executing Gaussian elimination with row-cyclic data distribution */
    double *x=gauss_cyclic(A, b, p, myrank);
	t2 = MPI_Wtime();

    /* Master rank processor write result x into output file and print the execution time */
    if (myrank==MASTER)
    {
    	output_matrix(fname_out, x, N);
    	printf("Wtime: %.9f \n", (t2-t1));
    }
    MPI_Finalize();
    return 0;
}
