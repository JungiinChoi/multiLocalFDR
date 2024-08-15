#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <headers.h>
#include <float.h>

void calcGradFloatC(double* gradA, double* gradB, double* influence, double* TermA, double* TermB, float* X, float* XW, float* grid, unsigned short int* YIdx, double* a, double* b, float gamma, float weight, float* delta, int N, int NIter, int M, int dim, int nH)
{
	float *grad_st_tmp = calloc(nH*(dim+1),sizeof(float));
	float *aGamma = malloc(dim*nH*sizeof(float)); 
	float *bGamma = malloc(nH*sizeof(float));
	int dimnH = dim*nH;
	int i,j,k;
	float factor = 1/gamma;
	float *gridLocal = malloc(dim*sizeof(float));	
	float TermALocal, TermBLocal;
	
	/* initialize some variables */
	for (k=0; k < dim; k++) {
		gridLocal[k] = grid[k];
	}

  	/* initialize some variables */
    for (i=0; i < nH; i++) {
        for (k=0; k < dim; k++) {
            aGamma[i*dim + k] = gamma*a[i*dim + k];
        }
        bGamma[i] = gamma*b[i];
		influence[i] = 0;
	}


	/* calculate gradient for samples */ 
	TermALocal = 0;	
	#pragma omp parallel
    {
   	  	float ftInnerMax;
        float sum_ft, sum_ft_inv;
        float *ft = calloc(nH,sizeof(float));
        float *ftInner = calloc(nH,sizeof(float));
        float *grad_ft_private = calloc(nH*(dim+1),sizeof(float));
		float Xtmp[dim];
		int *idxElements = malloc(nH*sizeof(int));
		int numElements, idxSave;
        float ftTmp;

        /* calculate gradient for samples X */
        #pragma omp for schedule(dynamic) private(i,k) reduction(+:TermALocal)
        for (j=0; j < NIter; j++) {
			for (k=0; k < dim; k++) {
				Xtmp[k] = X[j + k*N];
			}
			ftInnerMax = -FLT_MAX; numElements = 0;
            for (i=0; i < nH; i++) {
                ftTmp = bGamma[i] + aGamma[i*dim]*Xtmp[0];
                for (k=1; k < dim; k++) {
                    ftTmp += aGamma[i*dim + k]*Xtmp[k];
                }
                if (ftTmp > ftInnerMax - 25) {
                    if (ftTmp > ftInnerMax) {
                        ftInnerMax = ftTmp;
                    }
                    ftInner[numElements] = ftTmp;
                    idxElements[numElements++] = i;
                }
            }

            sum_ft = 0;
            /* calculate ft only for those entries that will be non-zero */
            for (i=0; i < numElements; i++) {
                ft[i] = expf(ftInner[i]-ftInnerMax);
                sum_ft += ft[i];
            }

            TermALocal += XW[j]*(ftInnerMax + logf(sum_ft))*factor;
//   			printf("%d: %.3f\n",j,XW[j]*(ftInnerMax + log(sum_ft))*factor);
   		    sum_ft_inv = 1/sum_ft*XW[j];

            /* update the gradient */
            for (i=0; i < numElements; i++) {
                idxSave = idxElements[i];
                for (k=0; k < dim; k++) {
                    grad_ft_private[idxSave + (k*nH)] += ft[i]*X[j+(k*N)]*sum_ft_inv;
                }
                grad_ft_private[idxSave + (dim*nH)] += ft[i]*sum_ft_inv;
            }
        }
		#pragma omp critical
        {
			for (i=0; i < nH*(dim+1); i++) {
            	gradA[i] += grad_ft_private[i];
            }
        }
		free(ft); free(ftInner); free(grad_ft_private); free(idxElements);
	}
	*TermA += TermALocal;

	/* Calculate gradient for grid points */ 
	TermBLocal = 0; 
	#pragma omp parallel 
	{
		float *Ytmp = calloc(dim,sizeof(float));
		float stInnerMax; float stInnerCorrection = 0;
		float sum_st, sum_st_inv, tmpVal, sum_st_inv2;
		float *st = calloc(nH,sizeof(float));
		float *stInner = calloc(nH,sizeof(float));
		float *grad_st_private = calloc(nH*(dim+1),sizeof(float));
		float *influencePrivate = calloc(nH,sizeof(float));
		int *idxElements = malloc(nH*sizeof(int));
		int numElements, idxSave;
		float stTmp;
		/* calculate gradient for grid points */
		#pragma omp for schedule(dynamic) private(i,k) reduction(+:TermBLocal)
		for (j=0; j < M; j++) {
			stInnerMax = -FLT_MAX; numElements = 0;
			for (k=0; k < dim; k++) {
				Ytmp[k] = gridLocal[k]+delta[k]*YIdx[k + j*dim];
			}
            for (i=0; i < nH; i++) {
                stTmp = bGamma[i] + aGamma[i*dim]*Ytmp[0];
                for (k=1; k < dim; k++) {
                    stTmp += (aGamma[i*dim + k]*Ytmp[k]);
                }
                if (stTmp > stInnerMax - 25) {
                    if (stTmp > stInnerMax) {
                        stInnerMax = stTmp;
                    }
                    stInner[numElements] = stTmp;
                    idxElements[numElements++] = i;
                }
            }

            sum_st = 0;
            /* calculate st only for those entries that wont be zero: exp^(-large number) */
            for (i=0; i < numElements; i++) {
                st[i] = expf(stInner[i]-stInnerMax);
                sum_st += st[i];
            }
			
			stInnerCorrection = expf(-stInnerMax*factor);
			tmpVal = powf(sum_st,-factor)*stInnerCorrection;
			
			TermBLocal += tmpVal; /*evalGrid[j] = tmpVal; */
			sum_st_inv2 = 1/sum_st;
			sum_st_inv = tmpVal*sum_st_inv2;
					

			for (i=0; i < numElements; i++) {
				idxSave = idxElements[i];
				influencePrivate[idxSave] += st[i]*sum_st_inv2;
		//		printf("%d: %.5f, %d, %d\n",j,st[i]*sum_st_inv2,numElements,idxElements[i]);
				st[i] *= sum_st_inv;
				grad_st_private[idxSave] += Ytmp[0]*st[i];
				for (k=1; k < dim; k++) {
					grad_st_private[idxSave + (k*nH)] += Ytmp[k]*st[i];
				}
				grad_st_private[idxSave + dimnH] += st[i];
			}
		}
		#pragma omp critical
		{
			for (i=0; i < nH; i++) {
				influence[i] += influencePrivate[i];
			}
			for (i=0; i < nH*(dim+1); i++) {
				grad_st_tmp[i] += grad_st_private[i];
			}
		}
		free(Ytmp); free(st); free(stInner); free(grad_st_private); free(influencePrivate); free(idxElements);
	} /* end of pragma parallel */

	*TermB += TermBLocal*weight;
	for (i=0; i < nH*(dim+1); i++) {
		gradB[i] -= grad_st_tmp[i]*weight;
	}
	
	free(grad_st_tmp); free(gridLocal); free(aGamma); free(bGamma);
}

