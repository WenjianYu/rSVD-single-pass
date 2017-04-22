#include "mkl_scalapack.h"
#include "matrix_funs_intel_mkl.h"

double get_seconds_frac(struct timeval start_timeval, struct timeval end_timeval){
    long secs_used, micros_used;
    secs_used=(end_timeval.tv_sec - start_timeval.tv_sec);
    micros_used= ((secs_used*1000000) + end_timeval.tv_usec) - (start_timeval.tv_usec);
    return (micros_used/1e6); 
}

/* initialize new matrix and set all entries to zero */
mat * matrix_new(int nrows, int ncols)
{
    mat *M = malloc(sizeof(mat));
    //M->d = (double*)mkl_calloc(nrows*ncols, sizeof(double), 64);
    M->d = (double*)calloc(nrows*ncols, sizeof(double));
    M->nrows = nrows;
    M->ncols = ncols;
    return M;
}


/* initialize new vector and set all entries to zero */
vec * vector_new(int nrows)
{
    vec *v = malloc(sizeof(vec));
    //v->d = (double*)mkl_calloc(nrows,sizeof(double), 64);
    v->d = (double*)calloc(nrows,sizeof(double));
    v->nrows = nrows;
    return v;
}


void matrix_delete(mat *M)
{
    //mkl_free(M->d);
    free(M->d);
    free(M);
}


void vector_delete(vec *v)
{
    //mkl_free(v->d);
    free(v->d);
    free(v);
}


/* copy contents of mat S to D  */
void matrix_copy(mat *D, mat *S){
    int i;
    //#pragma omp parallel for
    #pragma omp parallel shared(D,S) private(i) 
    {
    #pragma omp for 
    for(i=0; i<((S->nrows)*(S->ncols)); i++){
        D->d[i] = S->d[i];
    }
    }
}

/* initialize a random matrix */
void initialize_random_matrix(mat *M){
    int i,m,n;
    double val;
    m = M->nrows;
    n = M->ncols;
    float a=0.0,sigma=1.0;
    int N = m*n;
    float *r;
    VSLStreamStatePtr stream;
    
    r = (float*)malloc(N*sizeof(float));
   
    vslNewStream( &stream, BRNG,  time(NULL) );
    //vslNewStream( &stream, BRNG,  SEED );

    vsRngGaussian( METHOD, stream, N, r, a, sigma );

    // read and set elements
    #pragma omp parallel shared(M,N,r) private(i,val) 
    {
    #pragma omp parallel for
    for(i=0; i<N; i++){
        val = r[i];
        M->d[i] = val;
    }
    }
    
    free(r);
}


/* initialize new matrix and set all entries to zero  for float*/

void matrix_matrix_mult_row(mat *A, mat* B, mat* C){
    double alpha, beta;
    alpha = 1.0; beta = 0.0;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A->nrows, B->ncols, A->ncols, alpha, A->d, A->ncols, B->d, B->ncols, beta, C->d, C->ncols);
}

void matrix_transpose_matrix_mult_row(mat *A, mat* B, mat* C){
    double alpha, beta;
    alpha = 1.0; beta = 0.0;
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, A->ncols, B->ncols, A->nrows, alpha, A->d, A->ncols, B->d, B->ncols, beta, C->d, C->ncols);
}

/* C = A*B, A is a file on hard disk */
void matrix_matrix_mult_disk(FILE *A, mat *B, mat *C, int row, int col, int l){
    int row_size = l;
    int read_row_size = row_size;
    double alpha, beta;
    int i, j; // count
    int m=row, n=col, k=B->ncols;
    
    alpha = 1.0;
    beta = 0.0;
    // printf("matrix_matrix_mult_disk is running\n");
    float *M_f = (float*)malloc(read_row_size*n*sizeof(float));
    double *M = (double*)malloc(read_row_size*n*sizeof(double));

    struct timeval start_timeval_1, end_timeval_1;
    struct timeval start_timeval_2, end_timeval_2;
    double sum = 0;
    double time_1, time_2;
    gettimeofday(&start_timeval_1, NULL);

    for (i = 0; i < m; i += row_size){
        if (row_size > (m - i))
            read_row_size = m - i;
        gettimeofday(&start_timeval_2, NULL);   //time_2
        fread(M_f, sizeof(float), n*read_row_size, A);
        #pragma omp parallel shared(M, M_f,n,read_row_size) private(j) 
        {
            #pragma omp parallel for
            for(j=0; j<n*read_row_size; j++){
                M[j] = M_f[j];          //leixing zhuanhuan
                
            }
        }
        gettimeofday(&end_timeval_2, NULL);
        sum += get_seconds_frac(start_timeval_2 ,end_timeval_2);
       
        /* 1*n , n*k  =  1*k , all m*k */
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, read_row_size, B->ncols, n, alpha, M, n, B->d, B->ncols, beta, C->d+i*k, C->ncols);
        
    }
    gettimeofday(&end_timeval_1, NULL);
    time_1 = get_seconds_frac(start_timeval_1 ,end_timeval_1);
    
    time_2 = sum;
    printf("Time for reading data file_(fread-time1): %g second\n",time_2);
    printf("Time for matrix_matrix_mult: %g second\n", time_1);


    free(M_f);
    free(M);    
}
/* C = A^T*B ; column major */
/* n*m , m*k  =  n*k , all n*k */
void matrix_transpose_matrix_mult_disk(FILE *A, mat *B, mat *C, int row, int col, int l){
    int row_size = l;
    int read_row_size = row_size;

    double alpha, beta;
    int i, j; // count
    int m=row, n=col, k=B->ncols;
    float *M_f = (float*)malloc(read_row_size*n*sizeof(float));
    double *M = (double*)malloc(read_row_size*n*sizeof(double));;
    // printf("matrix_transpose_matrix_mult_disk is running\n");
    alpha = 1.0;
    beta = 1.0;
    // innitial C=0
    #pragma omp parallel shared(C) private(i) 
        {
            #pragma omp parallel for
            for(i=0; i < (C->nrows*C->ncols); i++){
                C->d[i] = 0.0;          
            }
        }
    

    struct timeval start_timeval_1, end_timeval_1;
    struct timeval start_timeval_2, end_timeval_2;
    double sum = 0;
    double time_1, time_2;
    gettimeofday(&start_timeval_1, NULL);

    for (i = 0; i < m; i += row_size){
        if (row_size > (m-i) )
            read_row_size = m-i;
        gettimeofday(&start_timeval_2, NULL);   //time_2
        fread(M_f, sizeof(float), n*read_row_size, A);
        // cblas_dcopy(k, B->d+i*k, 1, g_row, 1); //g_row = g[i];

        #pragma omp parallel shared(M,M_f,n,read_row_size) private(j) 
        {
            #pragma omp parallel for
            for(j=0; j < n * read_row_size; j++){
                M[j] = M_f[j];          
            }
        }
        gettimeofday(&end_timeval_2, NULL);
        sum += get_seconds_frac(start_timeval_2 ,end_timeval_2);
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, B->ncols, read_row_size, alpha, M, n, B->d+i*k, B->ncols, beta, C->d, C->ncols);

    }

    gettimeofday(&end_timeval_1, NULL);
    time_1 = get_seconds_frac(start_timeval_1 ,end_timeval_1);
    
    time_2 = sum;
    printf("Time for reading data file_(fread-time2): %g second\n",time_2);
    printf("Time for matrix_transpose_matrix_mult: %g second\n", time_1);

    free(M_f);
    free(M);
}

/* Performs [Q,R] = qr(M,'0') compact QR factorization 
M is mxn ; Q is mxn ; R is min(m,n) x min(m,n) */ 
void compact_QR_factorization(mat *M, mat *Q, mat *R){
    int i,j,m,n,k;
    m = M->nrows; n = M->ncols;
    k = min(m,n);

    mat *R_full = matrix_new(m,n);
    matrix_copy(R_full,M);
    //vec *tau = vector_new(n);
    vec *tau = vector_new(k);
    // get R
    //printf("get R..\n");
    //LAPACKE_dgeqrf(CblasColMajor, m, n, R_full->d, n, tau->d);
    LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, R_full->nrows, R_full->ncols, R_full->d, R_full->ncols, tau->d);
    
    for(i=0; i<k; i++){
        for(j=0; j<k; j++){
            if(j>=i){
                matrix_set_element(R,i,j,matrix_get_element(R_full,i,j));
            }
        }
    }

    // get Q
    matrix_copy(Q,R_full); 
    //printf("dorgqr..\n");
    LAPACKE_dorgqr(LAPACK_ROW_MAJOR, Q->nrows, Q->ncols, min(Q->ncols,Q->nrows), Q->d, Q->ncols, tau->d);


    // clean up
    matrix_delete(R_full);
    vector_delete(tau);
}


/* orth (Q)*/
void QR_factorization_getQ_inplace(mat *Q){
  
    int i,j,m,n,k;
    m = Q->nrows; n = Q->ncols;
    k = min(m,n);
    vec *tau = vector_new(k);

	/* do QR */	//sometime core dump, bug of MKL
    //LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, m, n, Q->d, n, tau->d);

	/* do QRCP */	//more stable, but more expensive
	printf("Warning: use QRCP to replace QR! (see line 269 of matrix_funs_intel_mkl.c)\n");
    int *jpvt = (int *)malloc(sizeof(int)*n);
     LAPACKE_dgeqpf(LAPACK_ROW_MAJOR, m, n, Q->d, n, jpvt, tau->d);
    free(jpvt);

    LAPACKE_dorgqr(LAPACK_ROW_MAJOR, m, n, k, Q->d, n, tau->d);
    vector_delete(tau);
}
/* M(:,inds) = Mc */
void matrix_set_selected_columns(mat *M, int *inds, mat *Mc){
    int i;
    vec *col_vec; 
    #pragma omp parallel shared(M,Mc,inds) private(i,col_vec) 
    {
    #pragma omp parallel for
    for(i=0; i<(Mc->ncols); i++){
        col_vec = vector_new(M->nrows); 
        matrix_get_col(Mc,i,col_vec);
        matrix_set_col(M,inds[i],col_vec);
        vector_delete(col_vec);
    }
    }
}

/* M(inds,:) = Mr */
void matrix_set_selected_rows(mat *M, int *inds, mat *Mr){ //modify
    int i;
    vec *row_vec; 
    #pragma omp parallel shared(M,Mr,inds) private(i,row_vec) 
    {
    #pragma omp parallel for
    for(i=0; i<(Mr->nrows); i++){
        row_vec = vector_new(M->ncols); 
        matrix_get_row(Mr,i,row_vec);
        matrix_set_row(M,inds[i],row_vec);
        vector_delete(row_vec);
    }
    }
}

/* Mc = M(:,inds) */
void matrix_get_selected_columns(mat *M, int *inds, mat *Mc){ //modify
    int i;
    vec *col_vec; 
    #pragma omp parallel shared(M,Mc,inds) private(i,col_vec) 
    {
    #pragma omp parallel for
    for(i=0; i<(Mc->ncols); i++){
        col_vec = vector_new(M->nrows);
        matrix_get_col(M,inds[i],col_vec);
        matrix_set_col(Mc,i,col_vec);
        vector_delete(col_vec);
    }
    }
}

/* extract column of a matrix into a vector */
void matrix_get_col(mat *M, int j, vec *column_vec){//modify
    int i;
    // unclear
    #pragma omp parallel shared(column_vec,M,j) private(i) 
    {//unclear
    #pragma omp parallel for
    for(i=0; i<M->nrows; i++){ 
        vector_set_element(column_vec,i,matrix_get_element(M,i,j));
    }
    }
}
/* extract row i of a matrix into a vector */
void matrix_get_row(mat *M, int i, vec *row_vec){//modify
    int j;
    #pragma omp parallel shared(row_vec,M,i) private(j) 
    {
    #pragma omp parallel for
    for(j=0; j<M->ncols; j++){ 
        vector_set_element(row_vec,j,matrix_get_element(M,i,j));
    }
    }
}
/* set column of matrix to vector */
void matrix_set_col(mat *M, int j, vec *column_vec){   //modify
    int i;
    #pragma omp parallel shared(column_vec,M,j) private(i) 
    {
    #pragma omp for
    for(i=0; i<M->nrows; i++){
        matrix_set_element(M,i,j,vector_get_element(column_vec,i));
    }
    }
}

/* put vector row_vec as row i of a matrix */
void matrix_set_row(mat *M, int i, vec *row_vec){  //modify
    int j;
    #pragma omp parallel shared(row_vec,M,i) private(j) 
    {
    #pragma omp parallel for
    for(j=0; j<M->ncols; j++){ 
        matrix_set_element(M,i,j,vector_get_element(row_vec,j));
    }
    }
}

/* set vector element */
void vector_set_element(vec *v, int row_num, double val){ //modify
    v->d[row_num] = val;
}

/* set element in column major format */
void matrix_set_element(mat *M, int row_num, int col_num, double val){ //modify
    M->d[row_num*(M->ncols) + col_num] = val;
}
/* get element in column major format */
double matrix_get_element(mat *M, int row_num, int col_num){ //modify
    return M->d[row_num*(M->ncols) + col_num];
}
double vector_get_element(vec *v, int row_num){
    return v->d[row_num];
}
/*********************Lijian***********************/


/* C = A*B & D = A^T*C */
void matrix_union_matrix_mult_disk_mem(FILE *A, mat *B, mat *C, mat *D, int row, int col, int row_size){
    int read_row_size = row_size;

    double alpha, beta ,gama;
    int i, j; // count
    int m=row, n=col, k=B->ncols;
    float *M_f = (float*)malloc(read_row_size*n*sizeof(float));
    double *M = (double*)malloc(read_row_size*n*sizeof(double));
    // double *g_row= (double*)malloc(k*sizeof(double)); //C's row vector'
    //printf("matrix_union_matrix_mult_disk_mem is running\n");
    
    alpha = 1.0;
    beta = 0.0;
    gama = 1.0;
    struct timeval start_timeval_1, end_timeval_1;
    struct timeval start_timeval_2, end_timeval_2;
    double sum = 0;
    double time_1, time_2;
    gettimeofday(&start_timeval_1, NULL);   //time_1
    //  #pragma omp parallel shared(D) private(i) 
    //     {
    //         #pragma omp parallel for
    //         for(i=0; i < (D->nrows*D->ncols); i++){
    //             D->d[i] = 0.0;          
    //         }
    //     }

    for (i = 0; i < m; i += row_size)
    {
        if (row_size > (m-i) )
            read_row_size = m-i;
        gettimeofday(&start_timeval_2, NULL);   //time_2
        fread(M_f, sizeof(float), n*read_row_size, A);
        #pragma omp parallel shared(M,M_f,n,read_row_size) private(j) 
        {
            #pragma omp parallel for
            for(j=0; j < n*read_row_size; j++){
                M[j] = M_f[j];          
            }
        }
        gettimeofday(&end_timeval_2, NULL);
        sum += get_seconds_frac(start_timeval_2 ,end_timeval_2);
        // C = A*D
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, read_row_size, B->ncols, n, alpha, M, n, B->d, B->ncols, beta, C->d+i*k, C->ncols);
        // B = A^T*C exchange B & C
        // cblas_dcopy(k, C->d+i*k, 1, g_row, 1); //g_row = g[i];
        // cblas_dger(CblasRowMajor, D->nrows, D->ncols, alpha, M, 1, g_row, 1, D->d, D->ncols); //A := alpha*x*y'+ A,
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, C->ncols, read_row_size, alpha, M, n, C->d+i*k, C->ncols, gama, D->d, D->ncols);
    }
    
    gettimeofday(&end_timeval_1, NULL);
    time_1 = get_seconds_frac(start_timeval_1 ,end_timeval_1);
    
    time_2 = sum;
    
    printf("Time for reading data file_(fread-time): %g second\n",time_2);
    printf("Time for matrix_union_matrix_mult: %g second\n", time_1);
    //printf("matrix_union_mem is %d KB\n",  getCurrentRSS()/1024);
    
 

    free(M_f);
    free(M);
}

//input A B C 
//output D E 
// D=A*B E=A^T*C
void matrix_union_matrix_mult_disk_mem_2(FILE *A, mat *B, mat *C, mat *D, mat*E, int row, int col, int row_size){
    int read_row_size = row_size;

    double alpha, beta ,gama;
    int i, j; // count
    int m=row, n=col, k=B->ncols;
    float *M_f = (float*)malloc(read_row_size*n*sizeof(float));
    double *M = (double*)malloc(read_row_size*n*sizeof(double));
    // double *g_row= (double*)malloc(k*sizeof(double)); //C's row vector'
    //printf("matrix_union_matrix_mult_disk_mem_2 is running\n");
    
    //matrix_copy(D,B); //D=B

    // float *a_g_mult=(float*)malloc(n*k*sizeof(float)); // ai * gi , n*k
    alpha = 1.0;
    beta = 0.0;
    gama = 1.0;
    struct timeval start_timeval_1, end_timeval_1;
    struct timeval start_timeval_2, end_timeval_2;
    double sum = 0;
    double time_1, time_2;
    gettimeofday(&start_timeval_1, NULL);   //time_1
    //  #pragma omp parallel shared(E) private(i) 
    //     {
    //         #pragma omp parallel for
    //         for(i=0; i < (E->nrows*E->ncols); i++){
    //             E->d[i] = 0.0;          
    //         }
    //     }

    for (i = 0; i < m; i += row_size)
    {
        if (row_size > (m-i) )
            read_row_size = m-i;
        gettimeofday(&start_timeval_2, NULL);   //time_2
        fread(M_f, sizeof(float), n*read_row_size, A);
        #pragma omp parallel shared(M,M_f,n,read_row_size) private(j) 
        {
            #pragma omp parallel for
            for(j=0; j < n*read_row_size; j++){
                M[j] = M_f[j];          
            }
        }
        gettimeofday(&end_timeval_2, NULL);
        sum += get_seconds_frac(start_timeval_2 ,end_timeval_2);
        // C = A*D
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, read_row_size, B->ncols, n, alpha, M, n, B->d, B->ncols, beta, D->d+i*k, D->ncols);
        // B = A^T*C exchange B & C
        // cblas_dcopy(k, C->d+i*k, 1, g_row, 1); //g_row = g[i];
        // cblas_dger(CblasRowMajor, D->nrows, D->ncols, alpha, M, 1, g_row, 1, D->d, D->ncols); //A := alpha*x*y'+ A,
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, C->ncols, read_row_size, alpha, M, n, C->d+(i*C->ncols), C->ncols, gama, E->d, E->ncols);
    }
    
    gettimeofday(&end_timeval_1, NULL);
    time_1 = get_seconds_frac(start_timeval_1 ,end_timeval_1);
    
    time_2 = sum;
    
    printf("Time for reading data file_(fread-time): %g second\n",time_2);
    printf("Time for matrix_union_matrix_mult2: %g second\n", time_1);

    free(M_f);
    free(M);
}
/*  k*n = k*k k*k n*k  */
void svd_row_cut (mat *A, mat *U, vec *E, mat * V)
{
    int m = A->nrows;
    int n = A->ncols;
    int i, j;
    // mat *A_in = matrix_new(m,n);;
    
    // matrix_copy(A_in, A);
    // printf("dong tai sheng qing\n");
    // double *u = (double*)malloc(m*m*sizeof(double));
    double *vt = (double*)malloc(n*m*sizeof(double));
    
    // printf("svd is running\n");
    // LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'S', m, n, A->d, n, E->d, U->d, m, vt, n, superb);
    
    
    // LAPACKE_dgesdd( int matrix_layout, char jobz, lapack_int m, lapack_int n, double* a, lapack_int lda, double* s, double* u, lapack_int ldu, double* vt, lapack_int ldvt ); 
    LAPACKE_dgesdd(LAPACK_ROW_MAJOR,'S', m, n, A->d, n, E->d, U->d, m, vt, n);
    //printf("Complete Lapack svd\n\n");

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            V->d[j * m + i] = vt[i * n+ j] ; 
        }
    }
    
    
    // printf("svd_row_cut is over\n");

    // matrix_delete(A_in);
    free(vt);
}

/* D = M(:,inds)' */
void matrix_get_selected_columns_and_transpose(mat *M, int *inds, mat *Mc){
    int i;
    vec *col_vec; 
    #pragma omp parallel shared(M,Mc,inds) private(i,col_vec) 
    {
    #pragma omp parallel for
    for(i=0; i<(Mc->nrows); i++){
        col_vec = vector_new(M->nrows);
        matrix_get_col(M,inds[i],col_vec);
        matrix_set_row(Mc,i,col_vec);
        vector_delete(col_vec);
    }
    }
}
void linear_solve_UTxb(mat *A, mat *b) {
    LAPACKE_dtrtrs(LAPACK_ROW_MAJOR, 'U', 'T', 'N',  //unclear
        b->nrows,
        b->ncols, 
        A->d,
        A->ncols,
        b->d,
        b->ncols
    );
}


/* C = beta*C + alpha*A(1:Anrows, 1:Ancols)[T]*B(1:Bnrows, 1:Bncols)[T] */
void submatrix_submatrix_mult_with_ab(mat *A, mat *B, mat *C, 
        int Anrows, int Ancols, int Bnrows, int Bncols, int transa, int transb, double alpha, double beta) {

    int opAnrows, opAncols, opBnrows, opBncols;
    if (transa == CblasTrans) {
        opAnrows = Ancols;
        opAncols = Anrows;
    } else {
        opAnrows = Anrows;
        opAncols = Ancols;
    }
    
    if (transb == CblasTrans) {
        opBnrows = Bncols;
        opBncols = Bnrows;
    } else {
        opBnrows = Bnrows;
        opBncols = Bncols;
    }
    
    if (opAncols != opBnrows) {
        printf("error in submatrix_submatrix_mult()");
        exit(0);
    }
    
    cblas_dgemm(CblasRowMajor, transa, transb, 
        opAnrows, opBncols, // m, n, 
        opAncols, // k
        alpha, A->d, A->ncols, // lda // modify
        B->d, B->ncols, // ldb
        beta, C->d, C->ncols // ldc
    );
    
}
void submatrix_submatrix_mult(mat *A, mat *B, mat *C, int Anrows, int Ancols, int Bnrows, int Bncols, int transa, int transb) {
    double alpha, beta;
    alpha = 1.0; beta = 0.0;
    submatrix_submatrix_mult_with_ab(A, B, C, Anrows, Ancols, Bnrows, Bncols, transa, transb, alpha, beta);
}
