#include "matrix_funs_intel_mkl.h"


// Algorithm 1
void rSVDbasic(char *filename,int m, int n, int k, int p, mat **Q, mat **B) {
    int i, j;
    FILE *fid;
    *Q = matrix_new(m, k);
    *B = matrix_new(k, n);
    mat *R, *G, *Bt;
    R = matrix_new(n, k);
    Bt = matrix_new(n, k);
    G = *Q;
    initialize_random_matrix(R);
    fid = fopen(filename, "rb");
    matrix_matrix_mult_disk(fid, R, G, m, n, k);
    fclose(fid);
    
    QR_factorization_getQ_inplace(G);
    fid = fopen(filename, "rb");
    matrix_transpose_matrix_mult_disk(fid, *Q, Bt, m, n, k);
    fclose(fid);
    for(i=0; i<(Bt->nrows); i++){
        for(j=0; j<(Bt->ncols); j++){
            matrix_set_element(*B,j,i,matrix_get_element(Bt,i,j)); 
        }
    }
    matrix_delete(Bt);
    matrix_delete(R);
}



// Algorithm 4
// caculate the memory of program
void rSVDsp(char *filename,int m, int n, int k, int kstep, int p, mat **Q, mat **B){
    int i, j;
    int row_size = k;
    mat *R, *G, *H;
    FILE *fid;
    //printf("rSVDsp_disk_mem is running\n");
    R = matrix_new(n, k);
    G = matrix_new(m, k);
    H = matrix_new(n, k);
    initialize_random_matrix(R);

    fid = fopen(filename, "rb");
    matrix_union_matrix_mult_disk_mem(fid, R, G, H, m, n, row_size);
    fclose(fid);
    
    if(p>0){   //power iteration  
		mat *Hi;
       	Hi = matrix_new(n, k);
    	for(i = 0; i < p; i++){
        	matrix_copy(Hi, H);
    		QR_factorization_getQ_inplace(Hi);

        	fid = fopen(filename, "rb");
        	matrix_union_matrix_mult_disk_mem(fid, Hi, G, H, m, n, row_size);
        	fclose(fid);
    	}
        matrix_delete(Hi);
	}
    // line 7-8

    *Q = matrix_new(m, k);
    *B = matrix_new(k, n);

     // initialize indices
    int *indices = (int*)calloc(k, sizeof(int));
    for (i = 0; i < k; i++) {
        indices[i] = i;
    }

    mat *omg, *g, *ht, *q, *r, *_q, *_r, *__r, 
        *Bomg, *QBomg, *omgtBt,
        *ytQ, *Qtq;

    omg = matrix_new(n, kstep);
    g = matrix_new(m, kstep);
    ht = matrix_new(kstep, n);
    q = matrix_new(m, kstep);
    _q = matrix_new(m, kstep);
    r = matrix_new(kstep, kstep);
    _r = matrix_new(kstep, kstep);
    __r = matrix_new(kstep, kstep);
    QBomg = matrix_new(m, kstep);
    // printf("QB is running\n");
    // line 10
    for (i = 0; i+kstep <= k; i+=kstep){
        matrix_get_selected_columns(R, indices+i, omg); 
        matrix_get_selected_columns(G, indices+i, g);
        matrix_get_selected_columns_and_transpose(H, indices+i, ht);
     
        // printf("i = %d\n", i);
        if (i == 0) {
            compact_QR_factorization(g, q, r);
            // compact_QR_factorization(g, q, r);
 
            linear_solve_UTxb(r, ht); //unclear
            matrix_set_selected_columns(*Q, indices+i, q);
            matrix_set_selected_rows(*B, indices+i, ht);

            
            
        } else {
            Bomg = matrix_new(i, kstep);
            omgtBt = matrix_new(kstep, i);
            Qtq = Bomg;    
            ytQ = omgtBt;   
            // line 9:
            // Bomg
            // printf("Bomg is running\n");
            submatrix_submatrix_mult(*B, omg, Bomg,
                    i, (*B)->ncols, omg->nrows, omg->ncols, CblasNoTrans, CblasNoTrans);
            // QBomg
            // printf("g is running\n");
            submatrix_submatrix_mult_with_ab(*Q, Bomg, g,
                    (*Q)->nrows, i, Bomg->nrows, Bomg->ncols, CblasNoTrans, CblasNoTrans, -1.0, 1.0);
            mat *y = g; 
            
           
            


            // line 10:
            compact_QR_factorization(y, q, r);
            
            // line 11:
            // Qtq
             submatrix_submatrix_mult(*Q, q, Qtq, 
                    (*Q)->nrows, i, q->nrows, q->ncols, CblasTrans, CblasNoTrans);
            submatrix_submatrix_mult_with_ab(*Q, Qtq, q,
                    (*Q)->nrows, i, Qtq->nrows, Qtq->ncols, CblasNoTrans, CblasNoTrans, -1.0, 1.0);
                    
            // _q, _r
            compact_QR_factorization(q, _q, _r);
 
            // line 12:
            // __r
            matrix_matrix_mult_row(_r, r, __r); //modify
            
            // line 13:
            // ytQ
            submatrix_submatrix_mult(y, *Q, ytQ,
                    y->nrows, y->ncols, (*Q)->nrows, i, CblasTrans, CblasNoTrans);
            submatrix_submatrix_mult_with_ab(ytQ, *B, ht,
                    ytQ->nrows, ytQ->ncols, i, (*B)->ncols, CblasNoTrans, CblasNoTrans, -1.0, 1.0);
            
            // omgtBt//line 16
            submatrix_submatrix_mult(omg, *B, omgtBt,
                    omg->nrows, omg->ncols, i, (*B)->ncols, CblasTrans, CblasTrans);
            submatrix_submatrix_mult_with_ab(omgtBt, *B, ht,
                    omgtBt->nrows, omgtBt->ncols, i, (*B)->ncols, CblasNoTrans, CblasNoTrans, -1.0, 1.0);
            linear_solve_UTxb(__r, ht);//line 17
            mat *b = ht;
            

            // line 12 - 13:
            matrix_set_selected_columns(*Q, indices+i, _q);


            matrix_set_selected_rows(*B, indices+i, b);

            
            matrix_delete(Bomg);
            matrix_delete(omgtBt);
        }
    }

    printf("rSVDsp total memory usage: %d KB\n", getCurrentRSS()/1024);
    
    matrix_delete(R);
    matrix_delete(G);
    matrix_delete(H);
    matrix_delete(omg);
    matrix_delete(r);
    matrix_delete(_r);
    matrix_delete(__r);
    matrix_delete(q);
    matrix_delete(_q);
    matrix_delete(g);
    matrix_delete(ht);            
    matrix_delete(QBomg);
    free(indices);
}

// Algorithm 2
void rSVD_exSP(char *filename,int m, int n, int k, int p, mat **Q,mat **_Q, mat **B) {
    int i, j;
    FILE *fid;
    *Q = matrix_new(m, k);
    *_Q = matrix_new(n, k);
    *B = matrix_new(k, k);
    mat *R, *_R , *Y, *_Y;
    R = matrix_new(n, k);
    _R = matrix_new(m, k);
    Y = matrix_new(m, k);
    _Y = matrix_new(n, k);

    initialize_random_matrix(R);
    initialize_random_matrix(_R);
    
    // fid = fopen(filename, "rb");
    // matrix_matrix_mult_disk_mem(fid, R, Y, m, n, k);
    // fclose(fid);

    fid = fopen(filename, "rb");
    matrix_union_matrix_mult_disk_mem_2(fid, R, _R, Y, _Y, m, n, k);
    fclose(fid);
    // power iteration
    // for (i = 0; i < p; i++) {
    //     fid = fopen(filename, "rb");
    //     matrix_transpose_matrix_mult_disk(fid, G, R, m, n);
    //     fclose(fid);
    //     fid = fopen(filename, "rb");
    //     matrix_matrix_mult_disk(fid, R, G, m, n);
    //     fclose(fid);
    // }
    matrix_copy(*Q, Y);
    matrix_copy(*(_Q), _Y);
    


    QR_factorization_getQ_inplace(*Q);
    QR_factorization_getQ_inplace(*_Q);
    
    mat *H = matrix_new(k,k);

    matrix_transpose_matrix_mult_row(_R, *Q, H);
    matrix_transpose_matrix_mult_row(_Y, *_Q, *B);
    
    // LAPACKE_dgesv (int matrix_layout , lapack_int n , lapack_int nrhs , double *a , lapack_int lda , lapack_int * ipiv , double * b , lapack_int ldb );
    int *ipiv = (int*)malloc(k*sizeof(int));
    LAPACKE_dgesv (LAPACK_ROW_MAJOR, k,k, H->d, k, ipiv, (*B)->d, k );
    
    printf("rSVD_exSP total memory usage: %d KB\n", getCurrentRSS()/1024);
    free(ipiv);
    matrix_delete(R);
    matrix_delete(_R);
    matrix_delete(Y);
    matrix_delete(_Y);
    matrix_delete(H);
}

