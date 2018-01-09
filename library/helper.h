#include <stdlib.h>
#include <stdio.h>
#include "mmio.h"
#include <stdbool.h>
#include <ctype.h>
#include <unistd.h>

/**
* @description      this funtion read the matrix from file and return the matrix by a double type pointer otherwise NULL
* @param    fname   input file name
* @param    n       matrix size
*/
double *input_matrix(const char* fname, int *n);

/**
* @description      this funtion write the result x into the file and return true if succeeded otherwise false
* @param    fname   input file name
* @param    x       calculation result
* @param    n       matrix size
*/
bool output_matrix(const char* fname, const double* x, int n);

/**
* @description      this funtion is used for handling command line options
* @param    argc    the number of arguments
* @param    argv    arguments
* @param    fin     input file name
* @param    fout    output file name
* @param    p1      the number of rows
* @param    p2      the number of columns
*/
bool handle_opt(int argc, char *argv[], char **fin, char **fout, int *p1, int *p2);

