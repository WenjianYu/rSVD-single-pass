#pragma once

#include <stdio.h>
#include "mkl.h"
#include "mkl_lapacke.h"
#include "mkl_vsl.h"

#include <time.h>
#include <sys/time.h> // for clock_gettime()


#define SEED    777
#define BRNG    VSL_BRNG_MCG31
#define METHOD  VSL_RNG_METHOD_GAUSSIAN_ICDF

typedef int BOOL;
#define true 1
#define false 0

#define min(x,y) (((x) < (y)) ? (x) : (y))
#define max(x,y) (((x) > (y)) ? (x) : (y))


typedef struct {
    int nrows, ncols;
    double * d;
} mat;


typedef struct {
    int nrows;
    double * d;
} vec;

double get_seconds_frac(struct timeval start_timeval, struct timeval end_timeval);
/* initialize new matrix and set all entries to zero */
mat * matrix_new(int nrows, int ncols);

/* initialize new vector and set all entries to zero */
vec * vector_new(int nrows);

void matrix_delete(mat *M);

void vector_delete(vec *v);

void matrix_matrix_mult_row(mat *A, mat* B, mat* C);
/* C= A*B */
void matrix_matrix_mult_disk(FILE *A, mat *B, mat *C, int row, int col, int k);
/* C= A^T*B */
void matrix_transpose_matrix_mult_disk(FILE *A, mat *B, mat *C, int row, int col, int k);
/* C = A*B & B = A^T*C */
void matrix_union_matrix_mult_disk(FILE *A, mat *B, mat *C, mat *D, int row, int col);

/* Performs [Q,R] = qr(M,'0') compact QR factorization 
M is mxn ; Q is mxn ; R is min(m,n) x min(m,n) */ 
void compact_QR_factorization_disk(mat *M, mat *Q, mat *R);
/* orth (Q) */
void QR_factorization_getQ_inplace_disk(mat *Q);

/* extract column of a matrix into a vector */
void matrix_get_col(mat *M, int j, vec *column_vec);
void matrix_get_row(mat *M, int i, vec *row_vec);
/* set column of matrix to vector */
void matrix_set_col(mat *M, int i, vec *column_vec);
/* put vector row_vec as row i of a matrix */
void matrix_set_row(mat *M, int i, vec *row_vec);
/* M(:,inds) = Mc */
void matrix_set_selected_columns(mat *M, int *inds, mat *Mc);
/* Mc = M(:,inds) */
void matrix_get_selected_columns(mat *M, int *inds, mat *Mc);
/* set vector element */
void vector_set_element(vec *v, int row_num, double val);
/* get element in column major format */
double matrix_get_element(mat *M, int row_num, int col_num);
double vector_get_element(vec *v, int row_num);
/* set element in column major format */
void matrix_set_element(mat *M, int row_num, int col_num, double val);


/* D = M(:,inds)' */
void matrix_get_selected_columns_and_transpose(mat *M, int *inds, mat *Mc);
void submatrix_submatrix_mult_disk(mat *A, mat *B, mat *C, int Anrows, int Ancols, int Bnrows, int Bncols, int transa, int transb);
void submatrix_submatrix_mult_with_ab_disk(mat *A, mat *B, mat *C, int Anrows, int Ancols, 
    int Bnrows, int Bncols, int transa, int transb, double alpha, double beta);

