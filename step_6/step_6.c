#include <stdlib.h>
#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h> 
#include "mmio.h"
#include "helper.h"

int N;
#define MASTER 0 //set master rank to be rank 0

/**
* @description      this funtion is used for getting the row r with maximum value on column k locally, and return r if succeeded otherwise -1
* @param    a   	matrix A
* @param    k    	target column k
* @param    p1    	the number of processors, column side
* @param    p2    	the number of processors, row side
* @param    me    	the rank of current processor
*/
int max_col_loc(double *a, int k, int p1, int p2,int me){
	double max=0;
	int i=0;
	int r = me / p1; //p2=4 r0~3

	while( r + i*p2 < k){
		i++;
	}

	int pos=r+i*p2;
	// printf("me:%d\tk:%d\tpos:%d\n", me, k, pos);

	if (pos>=N)
	{
		return -1;
	}
 
 	r=-1;
	for ( i = pos; i < N; i+=p2)
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
* @description      this funtion returns true if processor t is in the same columen group of k
* @param    k    	column k (a_kk)
* @param    p1    	the number of processors, column side
* @param    p2    	the number of processors, row side
* @param    t    	column t
*/
bool isInCoK(int k, int p1, int p2, int t){
    return (k%p1 == t%p1)?true:false;
}

/**
* @description      this funtion returns true if processor t is in the same row group of r
* @param    r    	row r
* @param    p1    	the number of processors, column side
* @param    p2    	the number of processors, row side
* @param    t    	column t
*/
bool isInRor(int r, int p1, int p2, int t){
	return (r%p2 == t/p1)?true:false;
}

/**
* @description      this funtion returns the column group owning column k belongs to
* @param    k    	column k
* @param    p1    	the number of processors, column side
* @param    p2    	the number of processors, row side
*/
int *Co(int k, int p1, int p2){
    int st = k % p1;
    // int co[p2];
    int *co = (int *)malloc( p2 * sizeof(*co));
    int i;
    for ( i = 0; i < p2; ++i)
    {
        co[i] = st + i * p1; 
    }
    return co;
}

/**
* @description      this funtion returns the row group owning row r belongs to
* @param    r    	row r
* @param    p1    	the number of processors, column side
* @param    p2    	the number of processors, row side
*/
int *Ro(int r, int p1, int p2){
    int *ro = (int *) malloc(p1*sizeof(*ro));
    int st = (r%p2) * p1;
    int i;
    for ( i = 0; i < p1; ++i)
    {
        ro[i] = st + i;
    }
    return ro;
}

/**
* @description      this funtion returns the column group to which a processor q belongs
* @param    q    	processor q
* @param    p1    	the number of processors, column side
* @param    p2    	the number of processors, row side
*/
int *Cop(int q, int p1, int p2){
    int i=0;
    int st = q % p1;
    int *cop = (int *)malloc( p2 * sizeof(*cop));
    for ( i = 0; i < p2; ++i)
    {
        cop[i] = st + (i * p1);
    }
    return cop;
}

/**
* @description      this funtion returns the row group to which a processor q belongs
* @param    q    	processor q
* @param    p1    	the number of processors, column side
* @param    p2    	the number of processors, row side
*/
int *Rop(int q, int p1, int p2){
    int i=0;
    int *rop = (int *) malloc(p1*sizeof(*rop));
    int st = q / p1;
    for ( i = 0; i < p1; ++i)
    {
        rop[i] = st*p1 + i; 
    }
    return rop;
}

/**
* @description      this funtion returns true if me is in the row group of processor q
* @param    q    	processor q
* @param    p1    	the number of processors, column side		
*/
bool isInRop(int me, int q, int p1){
	return (me/p1)==(q/p1)?true:false;
}

/**
* @description      this funtion returns true if me is in the column group of processor q
* @param    q    	processor q
* @param    p1    	the number of processors, column side
*/
bool isInCop(int me, int q, int p1){
	return (me%p1)==(q%p1)?true:false;
}

/**
* @description      this funtion returns the first processor in group (column group of column k)
* @param    k    	column k
* @param    p1    	the number of processors, column side
* @param    p2    	the number of processors, row side
*/
int grp_leader(int k, int p1, int p2 ){
	int *co = Co( k, p1, p2);
	return *(co+k%p2);
}

/**
* @description      this funtion is used for exchanging row r and k in matrix A and b locally
* @param    a   	matrix A
* @param    b    	matrix b
* @param    r     	row r
* @param    k    	row k
* @param    me    	the rank of current processor
*/
void exchange_row_loc(double *a, double *b, int r, int k, int p1, int p2, int me){
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
* @description      this funtion is used for copying row k in matrix A and b into buf locally
* @param    a   	matrix A
* @param    b    	matrix b
* @param    k    	row k
* @param    buf     buffer
* @param    me    	the rank of current processor
* @param    p1    	the number of processors, column side
*/
void copy_row_loc(double *a, double *b, int k, double *buf, int me, int p1){
	int i;
	for ( i = 0; i < N; ++i)
	{
		*(buf+i)=*(a + k*N +i);
	}
	*(buf + N)=*(b + k);
}

/**
* @description      this funtion computes the communication partner for processor me, which is the processor in row group ro belonging to the same column group as me.
* @param    ro   	row group ro
* @param    me    	the rank of current processor
* @param    p1    	the number of processors, column side
*/
int compute_partner(int *ro, int me, int p1){
	return *(ro + me%p1);
}

/**
* @description		this function computes the number of elements row k
* @param    k    	row k
*/
int compute_size(int k){
	return N-k+1;
}

/**
* @description      this funtion is used for exchanging row r in matrix A and b with buffer buf locally
* @param    a   	matrix A
* @param    b    	matrix b
* @param    r     	row r
* @param    buf     buffer
* @param    k    	row k
* @param    me    	the rank of current processor
* @param    psz    	the number of elements of the pivot row
*/
void exchange_row_buf(double *a, double *b, int r, double *buf, int k, int me, int psz){
	double temp;
	int i=k;
	
	for ( ; i < N; ++i)
	{
		temp = *(a + r*N +i);
		*(a + r*N +i) = *(buf+i);
		*(buf+i) = temp;
	}
	temp = *(b+r);
	*(b+r) = *(buf+N);
	*(buf+N) = temp;
}

/**
* @description      this funtion is used for computing the elimination factors l_ik for all elements a_ik owned by the processor. 
* @param    a   	matrix A
* @param    b    	matrix b
* @param    k    	row k
* @param    buf     buffer
* @param    p1    	the number of processors, column side
* @param    p2    	the number of processors, row side
* @param    me    	the rank of current processor
*/
double *compute_elim_fact_loc(double *a, double *b, int k, double *buf, int p1, int p2, int me){
	int i=k+1;
	double *l = (double *)malloc( N * sizeof(*l));
	while(i%p2 != me/p1)
		i++;
	for (; i < N; i+= p2){
		*(l+i)=(*(a+i*N+k))/(*(buf+k));
	}
	return l;
}

/**
* @description      	this funtion is used for computing the matrix elements
* @param    a   		matrix A
* @param    b    		matrix b
* @param    k    		row k
* @param    elim_buf	buffer storing elimination factors
* @param    buf     	buffer
* @param    me    		the rank of current processor
* @param    p1    		the number of processors, column side
* @param    p2    		the number of processors, row side
*/
void compute_local_entries(double *a, double *b, int k, double *elim_buf, double *buf, int me, int p1, int p2){
	int j;
	int i=k+1;
	while(i%p2 != me/p1){
		i++;
	}
	for (; i < N; i+=p2){
		j=k;
		while(j%p1!=me%p1){
			j++;
		}
		for (; j < N; j+=p1)
		{
			*(a+i*N+j) = *(a+i*N+j) - (*(elim_buf+i)) * (*(buf+j));
		}
		*(b+i) = *(b+i) - *(elim_buf+i) * *(buf+N);
	}
}

/**
* @description      this funtion executes Gaussian elimination with column-cyclic data distribution 
* @param    a   	matrix A
* @param    b    	matrix b
* @param    p1    	the number of processors, column side
* @param    p2    	the number of processors, row side, p2 will be 1 for column-cyclic data distribution 
* @param    p    	the number of processors
* @param    me    	the rank of current processor
*/
double * gauss_double_cyclic(double *a, double *b, int p1, int p2, int p, int me){
	double *x, *buf, *elim_buf;
	double sum, sum_g;
	int i, j, k, r, q, ql, size, buf_size, elim_size, psz, tag=42;
	struct 
	{
		double val;
		int pvtline;
	}z,y;
	MPI_Status status;

	x = (double *) malloc(N * sizeof(double));
	buf = (double *) malloc((N+1) * sizeof(double));
	elim_buf = (double *) malloc((N+1) * sizeof(double));

	MPI_Group world_group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
   	
   	/* Forward elimination*/  
	for ( k = 0; k < N-1 ; ++k)
	{
		/* if me is in the same column group of k */
		if (isInCoK(k, p1, p2, me))
		{
			/* get the local pivot row */
			r = max_col_loc(a, k, p1, p2, me);

			/* put the maximum |a_rk| with rank of current processor into z */
			z.pvtline = r;
			if(r==-1)
				z.val=0.;
			else
				z.val = fabs(*(a + r*N + k)); //a_rk

			/* create communicator: comm by communication group: group */
			MPI_Comm comm;
			MPI_Group group;
			MPI_Group_incl(world_group, p2, Co(k,p1,p2), &group);
			MPI_Comm_create_group(MPI_COMM_WORLD, group, 0, &comm);

			int prime_rank = -1, prime_size = -1;
			if (MPI_COMM_NULL != comm) {
			        MPI_Comm_rank(comm, &prime_rank);
			        MPI_Comm_size(comm, &prime_size);
			}			

			y.val=0;
			y.pvtline=-1;

			MPI_Reduce(&z, &y, 1, MPI_DOUBLE_INT, MPI_MAXLOC, k%p2, comm);
			MPI_Group_free(&group);
			if (MPI_COMM_NULL != comm) {
		    	MPI_Comm_free(&comm);
		    }
		}

		/* the first processor in column group of k broadcasts the maximum |a_rk|*/
		MPI_Bcast(&y, 1, MPI_DOUBLE_INT, grp_leader(k, p1, p2), MPI_COMM_WORLD);
		r = y.pvtline;

		/* pivot row and row k are on the same processor */
		/* column-cyclic data distribution, p2=1, k%p2 = r%p2, and isInRor(k, p1, p2, me) will always be true*/
		if ((k%p2 == r%p2))//k,r in same Ro
		{
			if (isInRor(k, p1, p2, me))
			{
				if (r!=k)
				{
					exchange_row_loc(a, b, r, k, p1, p2, me);
				}
				copy_row_loc(a,b,k,buf, me, p1);
			}
		}
		else{/* this branch will be executed in column-cyclic data distribution */
			if (isInRor(k, p1, p2, me))
			{
				copy_row_loc(a,b,k,buf, me, p1);
				
				q = compute_partner(Ro(r, p1, p2), me, p1);

				psz = compute_size(k);

				MPI_Send(buf+k, psz, MPI_DOUBLE, q, tag, MPI_COMM_WORLD);
			}
			else if (isInRor(r, p1, p2, me))//executing processor contains a part of the pivot row
			{
				q = compute_partner(Ro(k, p1, p2),me, p1);				
				psz = compute_size(k);
				MPI_Recv(buf+k, psz, MPI_DOUBLE, q, tag, MPI_COMM_WORLD, &status);
				exchange_row_buf(a, b, r, buf, k, me, psz);
			} 

		}

		/*broadcast of pivot row*/
		for (q = 0; q < p; ++q)
		{
			if (isInRor(r, p1, p2, q) &&  isInCop(me, q, p1)){
				ql = r % p2;//rank(q,Cop(q))
				buf_size = compute_size(k);

				/* create communicator: comm1 by communication group: group1 */
				MPI_Comm comm1;
				MPI_Group group1;
				MPI_Group_incl(world_group, p2, Cop(q,p1,p2), &group1);
				MPI_Comm_create_group(MPI_COMM_WORLD, group1, 0, &comm1);

				MPI_Bcast(buf+k, buf_size, MPI_DOUBLE, ql, comm1);
				MPI_Group_free(&group1);
				if (MPI_COMM_NULL != comm1) {
			    	MPI_Comm_free(&comm1);
			    }
			}
		}

		/* this if branch will not be executed as p2=1 for column-cyclic data distribution */
		if(((k%p2) != (r%p2)) && isInRor(k, p1, p2, me)){

			for ( i = 0; i < N; ++i)
			{
				*(a + k*N +i) = *(buf+i);
			}
			*(b+k)=*(buf+N);
		}

		/* if me is in the same column group of k, calculate elimination factors*/
		if (isInCoK(k, p1, p2, me))
		{
			elim_buf = compute_elim_fact_loc(a, b, k, buf, p1, p2, me);
		}

		/* broadcast of elimination factors */
		for (q = 0; q < p; ++q)
		{
			if (isInCoK(k, p1, p2, q) &&  isInRop(me, q, p1)){
				ql = (q) % p1 ;//rank(q,Cop(q))
				elim_size = N;

				/* create communicator: comm2 by communication group: group2 */
				MPI_Comm comm2;
				MPI_Group group2;
				MPI_Group_incl(world_group, p1, Rop(q,p1,p2), &group2);
				MPI_Comm_create_group(MPI_COMM_WORLD, group2, 0, &comm2);

				int prime_rank = -1, prime_size = -1;
				if (MPI_COMM_NULL != comm2) {
				        MPI_Comm_rank(comm2, &prime_rank);
				        MPI_Comm_size(comm2, &prime_size);
				}

				MPI_Bcast(elim_buf, elim_size, MPI_DOUBLE, ql, comm2);
				MPI_Group_free(&group2);
				if (MPI_COMM_NULL != comm2) {
			    	MPI_Comm_free(&comm2);
			    }

			}
		}
		/* computing the matrix elements */
		compute_local_entries(a, b, k, elim_buf, buf, me, p1, p2);
	}

	/* Backward substitution */
	for (k=N-1; k>=0; k--)
	{
		/* this condition will always be true for column-cyclic data distribution */
		if (isInRor(k, p1, p2, me))
		{
			sum=0.;
			j=k+1;
			while(j%p1!=me%p1){
				j++;
			}
			for ( ; j < N; j+=p1)
			{
				sum += *(a+k*N+j) * *(x+j);
			}

			/* create communicator: comm3 by communication group: group3 */
			MPI_Comm comm3;
			MPI_Group group3;
			MPI_Group_incl(world_group, p1, Rop(me,p1,p2), &group3);
			MPI_Comm_create_group(MPI_COMM_WORLD, group3, 0, &comm3);

			MPI_Allreduce(&sum, &sum_g, 1, MPI_DOUBLE, MPI_SUM, comm3);
			MPI_Group_free(&group3);
			if (MPI_COMM_NULL != comm3) {
		    	MPI_Comm_free(&comm3);
		    }

		    /* calculate x_k */
		    if (k%p1 == me%p1)
		    {
		    	*(x+k) = 1 / *(a+ k*N + k) * (*(b+k) - sum_g);
		    }
		}
		/* broadcast x_k among all processors*/
		MPI_Bcast(x+k, 1, MPI_DOUBLE, (k%p2)*p1+(k%p1), MPI_COMM_WORLD);
	}
	return x;
}

int main(int argc, char *argv[])
{
    int i, j, cur_msg_len, p, myrank;
    double t1, t2;
    
    /* Initialize the MPI, and get the number of processors into p, get the rank of current processor to myrank*/
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int p1=1, p2=1;
    char *fname_in;
    char *fname_out;
    /* get command line options, set input and output file names and the value of p1 and p2 */
    handle_opt( argc, argv, &fname_in, &fname_out, &p1, &p2);

	/*if the distribution of processors is not specified, let it be column-cyclic districution */
    if ( p1==-1 || p2 ==-1 )
    {
    	p1 = 1;
    	p2 = p;
    	if(myrank==MASTER)printf("no specification, p1xp2: %dx%d\n", p1, p2);
    }

    /* In step 6, p2 should be p and p1 should be 1 for column-cyclic districution */
    if (p1!=1)
    {
    	if(myrank==MASTER)printf("p1 need to be set to 1\n");
    	return 1;
    }

    /* check whether p1*p2 equals to p */
    if (p != p1*p2)
    {
    	if(myrank==MASTER)printf("ERROR: p not equals to p1xp2\n p:%d, p1xp2: %dx%d\n", p, p1, p2);
    	return 1;
    }

    double *A;

    /* only master rank processor read the matrix to A and the size of matrix size to N, then broadcast N to other processors*/
    if (myrank==MASTER)
    {
    	A=input_matrix(fname_in, &N);
    }

	MPI_Bcast(&N, 1, MPI_INT, MASTER, MPI_COMM_WORLD); 
    MPI_Barrier(MPI_COMM_WORLD);

    /* all elements in matrix b are set to 1 */
    double *b;
	b = (double *)malloc(N*sizeof(double));
    for ( i = 0; i < N; ++i)
	{
		*(b+i)=1;
	}
	if (myrank!= MASTER)
	{
		A = (double *)malloc(N*N*sizeof(double));
	}
    
    /* Master rank processor broadcast matrix A to all the other processors */
	MPI_Bcast(A, N * N , MPI_DOUBLE, MASTER, MPI_COMM_WORLD);   
    MPI_Barrier(MPI_COMM_WORLD);


    t1 = MPI_Wtime();
    /* Executing Gaussian elimination with column-cyclic data distribution */
    /* the definition of p1, p2 are not same in main() and gauss_double_cyclic() */
    /* the definition in main() follows the other steps whilst that in gauss_double_cyclic() follows the definition in the book */
    /*p1 is 1 (in main()) so that it's with column-cyclic data distribution */
    double *x=gauss_double_cyclic(A, b, p2, p1, p, myrank);
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
