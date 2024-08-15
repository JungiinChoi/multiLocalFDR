#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <headers.h>
#include <float.h>

void calcGradFastC(int* numEntries, int* elementList, int* maxElement, int* idxEntries, double* grad, double* influence, double* TermA, double* TermB, double* X, double* XW, double* grid, unsigned short int* YIdx, double* a, double* b, double gamma, double weight, double* delta, int N, int M, int dim, int nH)
{
	double *grad_st_tmp = calloc(nH*(dim+1),sizeof(double));
   	double *aGamma = malloc(dim*nH*sizeof(double));
    double *bGamma = malloc(nH*sizeof(double));
	int dimnH = dim*nH;
	int i,j,k;
	double factor = 1/gamma;
	double *gridLocal = malloc(dim*sizeof(double));	
	double TermALocal = 0, TermBLocal = 0;
    double epsCalcExp = -25; /* this is the maximal difference between the maximum and any other hyperplane in log-space; increase this value for more accuracy */

	/* initialize some variables */
	for (k=0; k < dim; k++) {
		gridLocal[k] = grid[k];
	}

	for (i=0; i < nH; i++) {
		for (j=0; j < dim; j++) {
			aGamma[i+(j*nH)] = gamma*a[i+(j*nH)];
		}
		bGamma[i] = gamma*b[i];
		influence[i] = 0;
	}

	for (i=0; i < nH*(dim+1); i++) {
		grad[i] = 0;
	}

	
	/* calculate gradient for samples */
	*TermA = 0;
	#pragma omp parallel
    {
   	  	double ftInnerMax;
        double sum_ft, sum_ft_inv;
        double *ft = calloc(nH,sizeof(double));
        double *ftInner = calloc(nH,sizeof(double));
        double *grad_ft_private = calloc(nH*(dim+1),sizeof(double));
		int *idxElements = malloc(nH*sizeof(int));
		int numElements, numElementsExp, idxSave, basicIdx;
        /* calculate gradient for samples */
        #pragma omp for private(i,k) reduction(+:TermALocal)
		for (j=0; j < N; j++) {
			ftInnerMax = -DBL_MAX;
			basicIdx = numEntries[j];
            numElements = numEntries[j+1] - basicIdx;
            for (i=0; i < numElements; i++) {
                ftInner[i] = bGamma[elementList[i+basicIdx]] + aGamma[elementList[i+basicIdx]]*X[j];
            }
            for (k=1; k < dim; k++) {
                for (i=0; i < numElements; i++) {
                    ftInner[i] += aGamma[elementList[i+basicIdx]+k*nH]*X[j + (k*N)];
                }
            }

			for (i=0; i < numElements; i++) {
				 if (ftInner[i] > ftInnerMax) {
                    ftInnerMax = ftInner[i];
                }
			}

			sum_ft = 0; numElementsExp = 0;
			/* calculate ft only for those entries that will be non-zero */
			for (i=0; i < numElements; i++) {
  				if (ftInner[i] - ftInnerMax > epsCalcExp) {
                    ft[numElementsExp] = exp(ftInner[i]-ftInnerMax);
                    idxElements[numElementsExp] = i;
                    sum_ft += ft[numElementsExp++];
                }

				/*ft[i] = exp(ftInner[i]-ftInnerMax);
				sum_ft += ft[i];*/
			}
			TermALocal += XW[j]*(ftInnerMax + log(sum_ft))*factor;
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
	*TermB = 0;
	#pragma omp parallel
	{
		double *Ytmp = calloc(dim,sizeof(double));
		double stInnerMax; double stInnerCorrection = 0;
		double sum_st, sum_st_inv, tmpVal, sum_st_inv2;
		double *stInner = calloc(nH,sizeof(double));
		double *grad_st_private = calloc(nH*(dim+1),sizeof(double));
		double *influencePrivate = calloc(nH,sizeof(double));
  		int *idxElements = malloc(nH*sizeof(int));
		int numElements, numElementsExp, idxSave, basicIdx;
		/* calculate gradient for grid points */
		#pragma omp for private(i,k) reduction(+:TermBLocal)
		for (j=0; j < M; j++) {
            basicIdx = numEntries[j+N];
            numElements = numEntries[j+1+N] - basicIdx;
			stInnerMax = -DBL_MAX;

			for (k=0; k < dim; k++) {
				Ytmp[k] = gridLocal[k]+delta[k]*YIdx[k + idxEntries[j]*dim];
			}

    		for (i=0; i < numElements; i++) {
                stInner[i] = bGamma[elementList[i+basicIdx]] + aGamma[elementList[i+basicIdx]]*Ytmp[0];
            }
            for (k=1; k < dim; k++) {
                for (i=0; i < numElements; i++) {
                    stInner[i] = stInner[i] + (aGamma[elementList[i+basicIdx]+k*nH]*Ytmp[k]);
                }
            }

   			for (i=0; i < numElements; i++) {
				if (stInner[i] > stInnerMax) {
                    stInnerMax = stInner[i];
                }
			}

			sum_st = 0; numElementsExp = 0;
			/* calculate st only for those entries that wont be zero */
	      	for (i=0; i < numElements; i++) {
	        	if (stInner[i] - stInnerMax > epsCalcExp) {
					stInner[numElementsExp] = exp(stInner[i]-stInnerMax);
					idxElements[numElementsExp] = i;
					sum_st += stInner[numElementsExp++];
                }
            }
			stInnerCorrection = exp(-stInnerMax*factor);
			tmpVal = pow(sum_st,-factor)*stInnerCorrection;
			
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
				influence[i] += influencePrivate[i];
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
