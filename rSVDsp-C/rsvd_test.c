#include "matrix_funs_intel_mkl.h"

/* parameter setting region */
//input data file
char filename[] = "exp3_2E3_2E3.dat"; 
int m = 2000;   				// row number of matrix             
int n = 2000;   				// column number of matrix          

//test values of rank, oversampling, and block size
int L[10] = {20,30,30, -1};		// rank parameter + oversampling
int RANK[10] = {16, 20, 24, -1};	// rank parameter

int KSTEP[10] = {10,-1};		// block size.  -1 is an ending token of array

//Other settings
BOOL outputflag= false;			//output U, S, V to disk file?
int round = 3;					//the times of running for each setting
int p4rSVDsp= 0; 				//0 or 1, only valid for rSVDsp alg.
/* end parameter setting */

char expname[50];

// algorithm 1
void rSVDbasic_test()
{
    // int ntrials = 3;
    int i, j;
    int k, ki, p=0;
    int rank;
    FILE *fid_U, *fid_V, *fid_S;
    int timesi;
    
    printf("----Begin rSVDbasic test!\n");
    for (ki = 0; L[ki]>0; ki++) 
    {   
        k=L[ki]; 
        rank = RANK[ki];

        printf("k = %d\n", rank);
                
        struct timeval start_timeval, end_timeval;
        double sum = 0;
        double time;
        for (timesi = 0; timesi < round; timesi++)
        {
            printf("\nround No. %d\n", timesi);
            mat *Q, *B;
            mat *U, *V;
            mat *Utrue;
            vec *E;
            U = matrix_new(k,k);
            V = matrix_new(n,k);
            E = vector_new(min(k,n));
            Utrue = matrix_new(m, k);
            gettimeofday(&start_timeval, NULL);    
            rSVDbasic(filename, m, n, k, p, &Q, &B);
                 
            svd_row_cut (B, U, E, V);
            matrix_matrix_mult_row(Q, U, Utrue);
            gettimeofday(&end_timeval, NULL);
            if (outputflag && timesi == 0)
            {
                char U_buffer[100];
                char V_buffer[100];
                char S_buffer[100];
                sprintf(U_buffer,"%s_rSVDbasic_k=%d_U.dat", expname,rank);
                sprintf(V_buffer,"%s_rSVDbasic_k=%d_V.dat", expname,rank);
                sprintf(S_buffer,"%s_rSVDbasic_k=%d_S.dat", expname,rank);
                // printf("%s",U_buffer);
                fid_U = fopen(U_buffer,"w");//mark
                fid_V = fopen(V_buffer,"w");//mark
                fid_S = fopen(S_buffer,"w");//mark

                double *u = (double*)malloc(rank*m*sizeof(double));
                double *v = (double*)malloc(n*rank*sizeof(double));
                        
                #pragma omp parallel shared(rank,m, u, Utrue) private(i,j) 
                {
                    #pragma omp parallel for
                    for (i = 0; i < m; i++)
                    {
                        for (j = 0; j< rank; j++)
                        {
                            u[i*rank+j] = Utrue->d[i*Utrue->ncols + j];
                        }
                    }
                }
                #pragma omp parallel shared(rank,n, v, V) private(i,j) 
                {
                    #pragma omp parallel for
                    for (i = 0; i< n; i++)
                    {
                        for (j = 0; j<rank;j++)
                        {
                            v[i*rank+j] = V->d[i*V->ncols +j];
                        }
                    }
                }
                        
                fwrite(u, sizeof(double), m*rank, fid_U);
                fwrite(v, sizeof(double), n*rank, fid_V);
                fwrite(E->d, sizeof(double), E->nrows, fid_S);
                fclose(fid_U);
                fclose(fid_V);
                fclose(fid_S);
                free(u);
                free(v);
            }
                    
            matrix_delete(U);
            matrix_delete(V);
            vector_delete(E);
            matrix_delete(Q); 
            matrix_delete(B);
            matrix_delete(Utrue);
            sum += get_seconds_frac(start_timeval,end_timeval);
        }
                
        time = sum / round;
        printf("\nAverage runtime of rSVDbasic alg.: %f second\n\n", time);
    }
}


// algorithm 4
void rSVDsp_test()
{
    // int ntrials = 3;
    int i, j;
    int  p, k, kstep, ni, ki, kstepi;
    int rank;
    FILE *fid_U, *fid_V, *fid_S;
    int timesi;
    
    p= p4rSVDsp;	//0 default, 1 with an extra pass.

    printf("----Begin rSVDsp test!\n");
    for (ki = 0; L[ki]>0; ki++) 
    {   
        k=L[ki]; 
        rank = RANK[ki];

        for (kstepi=0; KSTEP[kstepi]>0; kstepi++) 
        {   
            kstep= KSTEP[kstepi];
            printf("k = %d, kstep = %d, p=%d\n", rank, kstep, p);
                
            struct timeval start_timeval, end_timeval;
            double sum = 0;
            double time;
            for (timesi = 0; timesi < round; timesi++)
            {
                printf("\nround No. %d\n", timesi);
                mat *Q, *B;
                mat *U, *V;
                mat *Utrue;
                vec *E;
                U = matrix_new(k,k);
                V = matrix_new(n,k);
                E = vector_new(min(k,n));
                Utrue = matrix_new(m, k);
                gettimeofday(&start_timeval, NULL);    
                rSVDsp(filename,m, n, k, kstep, p, &Q, &B);
                svd_row_cut (B, U, E,  V);
                matrix_matrix_mult_row(Q, U, Utrue);
                gettimeofday(&end_timeval, NULL);
                if (outputflag && timesi == 0)
                {
                       
                    char U_buffer[100];
                    char V_buffer[100];
                    char S_buffer[100];
                    sprintf(U_buffer,"%s_rSVDsp_k=%d_U.dat", expname,rank);
                    sprintf(V_buffer,"%s_rSVDsp_k=%d_V.dat", expname,rank);
                    sprintf(S_buffer,"%s_rSVDsp_k=%d_S.dat", expname,rank);

                    fid_U = fopen(U_buffer,"w");
                    fid_V = fopen(V_buffer,"w");
                    fid_S = fopen(S_buffer,"w");
                    double *u = (double*)malloc(rank*m*sizeof(double));
                    double *v = (double*)malloc(n*rank*sizeof(double));
                        
                    #pragma omp parallel shared(rank,m, u, Utrue) private(i,j) 
                    {
                        #pragma omp parallel for
                        for (i = 0; i < m; i++)
                        {
                            for (j = 0; j< rank; j++)
                            {
                                u[i*rank+j] = Utrue->d[i*Utrue->ncols + j];
                            }
                        }
                    }
                    #pragma omp parallel shared(rank,n, v, V) private(i,j) 
                    {
                        #pragma omp parallel for
                        for (i = 0; i< n; i++)
                        {
                            for (j = 0; j<rank;j++)
                            {
                                v[i*rank+j] = V->d[i*V->ncols +j];
                            }
                        }
                    }
                    //printf("PCA_svd_mem takes up %d\n", getCurrentRSS()/1024);
                    fwrite(u, sizeof(double), m*rank, fid_U);
                    fwrite(v, sizeof(double), n*rank, fid_V);
                    fwrite(E->d, sizeof(double), E->nrows, fid_S);
                    fclose(fid_U);
                    fclose(fid_V);
                    fclose(fid_S);
                    free(u);
                    free(v);
                 }
                // fid_time = fopen("result_mem/result1.txt","a");
                //printf("PCA_mem takes up %d\n", getCurrentRSS()/1024);
                // fclose(fid_time);
                matrix_delete(U);
                matrix_delete(V);
                vector_delete(E);
                matrix_delete(Q); 
                matrix_delete(B);
                matrix_delete(Utrue);
                sum += get_seconds_frac(start_timeval,end_timeval);;
            }
                
            time = sum / round;
            printf("\nAverage runtime of rSVDsp alg.: %f second\n\n", time);
        }
    }
}

//algorithm 2
void rSVD_exSP_test()
{
    // int ntrials = 3;
    int i, j;
    int  p=2, k, kstep, ni, ki;
    int rank;
    FILE *fid_U, *fid_V, *fid_S;
    int timesi;
    
    printf("----Begin rSVD_exSP test!\n");

    for (ki = 0; L[ki]>0; ki++) 
    {   
        k=L[ki]; 
        rank = RANK[ki];
        // fid_time = fopen("result_mem/result2.txt","a"); //mark
        // fprintf(fid_time,"rank = %d, k = %d\n", rank, k);
        //printf("rank = %d, l = %d\n", rank, k);
        // fclose(fid_time);

        printf("k = %d\n", rank);
                
        struct timeval start_timeval, end_timeval;
        double sum = 0;
        double time;
        for (timesi = 0; timesi < round; timesi++)
        {
            printf("\nround No. %d\n", timesi);
            mat *Q, *B, *_Q;
            mat *U, *V;
            mat *Utrue;
            mat *H;
            vec *E;
            U = matrix_new(k,k);
            V = matrix_new(k,k);
            E = vector_new(k);
            Utrue = matrix_new(m, k);
            H = matrix_new(n, k);
            gettimeofday(&start_timeval, NULL);    
            rSVD_exSP(filename,m, n, k, p, &Q, &_Q, &B);

            svd_row_cut (B, U, E, V);
            matrix_matrix_mult_row(Q, U, Utrue);
            matrix_matrix_mult_row(_Q, V, H);
            gettimeofday(&end_timeval, NULL);
            if (outputflag && timesi == 0)
            {
                char U_buffer[100];
                char V_buffer[100];
                char S_buffer[100];
                sprintf(U_buffer,"%s_rSVD_exSP_k=%d_U.dat", expname,rank);
                sprintf(V_buffer,"%s_rSVD_exSP_k=%d_V.dat", expname,rank);
                sprintf(S_buffer,"%s_rSVD_exSP_k=%d_S.dat", expname,rank);

                fid_U = fopen(U_buffer,"w");//mark
                fid_V = fopen(V_buffer,"w");//mark
                fid_S = fopen(S_buffer,"w");//mark
                double *u = (double*)malloc(rank*m*sizeof(double));
                double *v = (double*)malloc(n*rank*sizeof(double));
                
                #pragma omp parallel shared(rank,m, u, Utrue) private(i,j) 
                {
                    #pragma omp parallel for
                    for (i = 0; i < m; i++)
                    {
                        for (j = 0; j< rank; j++)
                        {
                            u[i*rank+j] = Utrue->d[i*Utrue->ncols + j];
                        }
                    }
                }
                #pragma omp parallel shared(rank,n, v, H) private(i,j) 
                {
                    #pragma omp parallel for
                    for (i = 0; i< n; i++)
                    {
                        for (j = 0; j<rank;j++)
                        {
                            v[i*rank+j] = H->d[i*H->ncols +j];
                        }
                    }
                }
                //printf("PCA_svd_mem takes up %d\n", getCurrentRSS()/1024);
                fwrite(u, sizeof(double), m*rank, fid_U);
                fwrite(v, sizeof(double), n*rank, fid_V);
                fwrite(E->d, sizeof(double), E->nrows, fid_S);
                fclose(fid_U);
                fclose(fid_V);
                fclose(fid_S);
                free(u);
                free(v);
             }
            // fid_time = fopen("result_mem/result1.txt","a");
            //printf("takes up %d\n", getCurrentRSS()/1024);
            // fclose(fid_time);
            matrix_delete(U);
            matrix_delete(V);
            vector_delete(E);
            matrix_delete(Q); 
            matrix_delete(B);
            matrix_delete(_Q);
            matrix_delete(Utrue);
            matrix_delete(H);
                    
            sum += get_seconds_frac(start_timeval,end_timeval);;
            //sum = get_seconds_frac(start_timeval,end_timeval);
            //printf("%f\n", get_seconds_frac(start_timeval,end_timeval));
        }
        time = sum / round;
        printf("\nAverage runtime of rSVD_exSP alg.: %f second\n\n", time);
    }
}



int main() {
	char * ptr1, *ptr2;
	int len, i;
	ptr1= strrchr(filename, '/');
	ptr2= strchr(filename, '.');
	if(ptr1== NULL){
	 	i=0;
		while(filename[i]!='.'){
			expname[i]= filename[i];
			i++;
		}
		expname[i]='\0';
	}
	else{
		len= ptr2-ptr1-1;
		for(i=0; i<len; i++){
			expname[i]= ptr1[i+1];
		}
		expname[i]='\0';
	}

	printf("***************************************************************\n");
	printf("Test single-pass SVD algorithm (rSVDsp) and other counterparts.\n");
	printf("***************************************************************\n\n");
	printf("Input matrix file: %s\n\n", filename);

    rSVDbasic_test();
    rSVD_exSP_test();
    rSVDsp_test();
    
}
