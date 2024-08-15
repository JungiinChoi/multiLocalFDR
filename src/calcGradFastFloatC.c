#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <headers.h>
#include <float.h>

#define epsCalcExp -25

void calcGradFastFloatC(int* numEntries, int* elementList, int* maxElement, int* idxEntries, double* grad, double* influence, double* TermA, double* TermB, float* X, float* XW, float* grid, unsigned short int* YIdx, double* a, double* b, float gamma, float weight, float* delta, int N, int M, int dim, int nH)
{
	float *grad_st_tmp = calloc(nH*(dim+1),sizeof(float));
   	float *aGamma = malloc(dim*nH*sizeof(float));
    float *bGamma = malloc(nH*sizeof(float));
	int dimnH = dim*nH;
	int i,j,k;
	float factor = 1/gamma;
	float *gridLocal = malloc(dim*sizeof(float));	
	float TermALocal = 0, TermBLocal = 0;

	/* initialize some variables */
	for (k=0; k < dim; k++) {
		gridLocal[k] = grid[k];
	}

	for (i=0; i < nH; i++) {
		for (k=0; k < dim; k++) {
			aGamma[i*dim + k] = gamma*a[i*dim + k];
		}
		bGamma[i] = gamma*b[i];
		influence[i] = 0;
	}

	for (i=0; i < nH*(dim+1); i++) {
		grad[i] = 0;
	}

	
	/* calculate gradient for samples */
	#pragma omp parallel
    {
   	  	float ftInnerMax, ftTmp;
        float sum_ft, sum_ft_inv;
        float *ft = calloc(nH,sizeof(float));
        float *ftInner = calloc(nH,sizeof(float));
        float *grad_ft_private = calloc(nH*(dim+1),sizeof(float));
		int *idxElements = malloc(nH*sizeof(int));
		int numElements, numElementsExp, idxSave, basicIdx;
        /* calculate gradient for samples */
        #pragma omp for private(i,k) reduction(+:TermALocal) 
		for (j=0; j < N; j++) {
			ftInnerMax = -FLT_MAX;
			basicIdx = numEntries[j];
            numElements = numEntries[j+1] - basicIdx;
   			sum_ft = 0; numElementsExp = 0;
         	for (i=0; i < numElements; i++) {
            	ftTmp = bGamma[elementList[i+basicIdx]] + aGamma[elementList[i+basicIdx]*dim]*X[j];
	           	for (k=1; k < dim; k++) {
                	ftTmp += aGamma[elementList[i+basicIdx]*dim+k]*X[j + (k*N)];
            	}
				if (ftTmp - ftInnerMax > epsCalcExp) {
					if (ftTmp > ftInnerMax) {
						ftInnerMax = ftTmp;
					}
					ft[numElementsExp] = ftTmp;
		         	idxElements[numElementsExp++] = i;
				}
			}

			/* calculate ft only for those entries that will be non-zero */
			for (i=0; i < numElementsExp; i++) {
                ft[i] = expf(ft[i]-ftInnerMax);
                sum_ft += ft[i];
			}

			TermALocal += XW[j]*(ftInnerMax + logf(sum_ft))*factor;
			sum_ft_inv = 1/sum_ft;
		
			for (i=0; i < numElementsExp; i++) {
		        idxSave = elementList[idxElements[i]+basicIdx];
				for (k=0; k < dim; k++) {
					grad_ft_private[idxSave + (k*nH)] += XW[j]*ft[i]*X[j+(k*N)]*sum_ft_inv;
				}
				grad_ft_private[idxSave + (dim*nH)] += XW[j]*ft[i]*sum_ft_inv;
			}
		}
    	#pragma omp critical
        {
			for (i=0; i < nH*(dim+1); i++) {
                grad[i] += grad_ft_private[i];
            }
        }
		free(ft); free(ftInner); free(grad_ft_private); free(idxElements);
	}
	*TermA = TermALocal;

	/* Calculate gradient for grid points */
	#pragma omp parallel
	{
		float *Ytmp = calloc(dim,sizeof(float));
		float stInnerMax; float stInnerCorrection = 0;
		float sum_st, sum_st_inv, tmpVal, sum_st_inv2;
		float *stInner = calloc(nH,sizeof(float));
		float stTmp;
		float *grad_st_private = calloc(nH*(dim+1),sizeof(float));
		float *influencePrivate = calloc(nH,sizeof(float));
  		int *idxElements = malloc(nH*sizeof(int));
		int numElements, numElementsExp, idxSave, basicIdx;
		/* calculate gradient for grid points */
		#pragma omp for private(i,k) reduction(+:TermBLocal)
		for (j=0; j < M; j++) {
            basicIdx = numEntries[j+N];
            numElements = numEntries[j+1+N] - basicIdx;
			stInnerMax = -FLT_MAX;
	
			for (k=0; k < dim; k++) {
				Ytmp[k] = gridLocal[k]+delta[k]*YIdx[k + idxEntries[j]*dim];
			}

   			sum_st = 0; numElementsExp = 0;
 			for (i=0; i < numElements; i++) {
                stTmp = bGamma[elementList[i+basicIdx]] + aGamma[elementList[i+basicIdx]*dim]*Ytmp[0];
				for (k=1; k < dim; k++) {
					stTmp += aGamma[elementList[i+basicIdx]*dim+k]*Ytmp[k];
				}
				if (stTmp > stInnerMax) {
					stInnerMax = stTmp;
				}
				stInner[i] = stTmp;
			}

			/* calculate st only for those entries that wont be zero */
	      	for (i=0; i < numElements; i++) {
	        	if (stInner[i] - stInnerMax > epsCalcExp) {
					stInner[numElementsExp] = expf(stInner[i]-stInnerMax);
					idxElements[numElementsExp] = i;
					sum_st += stInner[numElementsExp++];
                }
            }
			stInnerCorrection = expf(-stInnerMax*factor);
			tmpVal = powf(sum_st,-factor)*stInnerCorrection;
			
			TermBLocal += tmpVal; //evalGrid[idxEntries[j]] = tmpVal;
			sum_st_inv2 = 1/sum_st;
			sum_st_inv = tmpVal*sum_st_inv2;
	
			for (i=0; i < numElementsExp; i++) {
			    idxSave = elementList[idxElements[i]+basicIdx];
				influencePrivate[idxSave] += stInner[i]*sum_st_inv2;
				stInner[i] *= sum_st_inv;
				grad_st_private[idxSave] += Ytmp[0]*stInner[i];
				for (k=1; k < dim; k++) {
					grad_st_private[idxSave + (k*nH)] += Ytmp[k]*stInner[i];
				}
				grad_st_private[idxSave + dimnH] += stInner[i];
			}
		}
		#pragma omp critical
		{
			for (i=0; i < nH; i++) {
				influence[i] += (double) influencePrivate[i];
			}
			for (i=0; i < nH*(dim+1); i++) {
				grad_st_tmp[i] += grad_st_private[i];
			}
		}
		free(Ytmp); free(stInner); free(grad_st_private); free(influencePrivate); free(idxElements);
	} /* end of pragma parallel */
	*TermB = TermBLocal*weight;
	for (i=0; i < nH*(dim+1); i++) {
		grad[i] -= grad_st_tmp[i]*weight;
	}
	free(grad_st_tmp); free(gridLocal); free(aGamma); free(bGamma); 
}
