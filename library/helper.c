#include "helper.h"

double *input_matrix(const char* fname, int *n){
    FILE *f;
    MM_typecode matcode;
    int M, N;
    int i,j;
    double *a=NULL;
    double *matrix;

    /* Check if the file is successfully opened */
    if ((f = fopen(fname, "r")) == NULL){
        printf("e1\n");
        return NULL;
    }
 
    /* read and check Matrix Market banner */
    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("mm_read_unsymetric: Could not process Matrix Market banner ");
        printf(" in file [%s]\n", fname);
        return NULL;
    }

    /* check whether the matrix is complex, sparse matrix */
    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && mm_is_sparse(matcode) )
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        return NULL;
    }

    /* get the matrix size the saved in M, N */
    if (mm_read_mtx_array_size(f, &M, &N) != 0) {
        printf("Could not read matrix's size.\n");
        return NULL;
    }

    /* Check if the matrix is square matrix */
    if(M!=N){
        printf("M != N\n");
        return NULL;
    }
    *n = M;

    /* Check if the matrix is dense, then read matrix and saved in matrix*/
    if (mm_is_dense(matcode))
    {
        matrix = (double *)malloc(*n * *n * sizeof(double));
        for ( j = 0; j < *n; ++j)
        {
            for ( i = 0; i < *n; ++i)
            {
                fscanf(f, "%lg", &matrix[i* *n+j]);
            }
        }
        return matrix;
    }
    printf("e2\n");
    return NULL;
}

bool output_matrix(const char* fname, const double* x, int n){
    FILE *stdout;
    if ((stdout = fopen(fname, "w")) == NULL){
        printf("e1\n");
        return NULL;
    }
    MM_typecode matcode;
    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_array(&matcode);
    mm_set_real(&matcode);

    mm_write_banner(stdout, matcode); 
    mm_write_mtx_array_size(stdout, n, 1);
    int i;
    for (i=0; i<n; i++)
        fprintf(stdout, "%.10lf\n", *(x+i));

    return true;
}

bool handle_opt(int argc, char *argv[], char **fin, char **fout, int *p1, int *p2){
    int opt, index;
    /* get options */
    opt = getopt( argc, argv, "i:o:g:" );
    char *p1p2=NULL;

    while( opt != -1 ) {
        switch(opt){
            case 'i':{
                *fin = optarg;
                break;
            }
            case 'o':{
                *fout = optarg;
                break;
            }
            case 'g':{
                p1p2 = optarg;
                break;
            }
            case '?':{
                if (optopt == 'i'){
                    fprintf (stderr, "Option -%c requires an argument.\n missing input file", optopt);
                }
                else if (opt == 'o')
                {
                    fprintf (stderr, "Option -%c requires an argument.\n missing output file", optopt);
                }
                else if (opt == 'g')
                {
                    fprintf (stderr, "Option -%c requires an argument.\n missing cores distribution", optopt);
                }

                else if (isprint (optopt))
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else
                    fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);

                return false;
            }
            default:{
                printf("switch default, opt:%c\n", (char) opt);
                return false;
            }
        }
        /*  get next option */
        opt = getopt( argc, argv, "i:o:g:" );
    }

    for (index = optind; index < argc; index++){
        printf("Non-option argument %s\n", argv[index]);
    }
        
    if (*fin==NULL)
    {
        *fin = "matrix_0008.mm";
        printf("no input file\n");
    }
    if (*fout==NULL)
    {
        *fout = "result2.mm";
        printf("no input file\n");
    }

    if (p1p2 == NULL) {
        *p1 = -1;
        *p2 = -1;
    }
    else{
        if (sscanf(p1p2, "%dx%d", p1, p2) != 2) {//ascll *
            printf("Could not parse -g option with argument: %s\n", p1p2);
            *p1 = -1;
            *p2 = -1;
        }
    }

    return true;
}




