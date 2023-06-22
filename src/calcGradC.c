#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <headers.h>
#include <float.h>
#include <limits.h>

void calcGradC(double* gradA, double* gradB, double* influence, double* TermA, double* TermB, double* X, double* XW, double* grid, unsigned short int* YIdx, int *numPointsPerBox, double* boxEvalPoints, unsigned short int *XToBox, int numBoxes, double* a, double* b, double gamma, double weight, double* delta, int N, int M, int dim, int nH, int MBox) {
	double *grad_st_tmp = calloc(nH*(dim+1),sizeof(double));
	double *aGamma = malloc(dim*nH*sizeof(double)); 
	double *bGamma = malloc(nH*sizeof(double));
	int dimnH = dim*nH;
	int i,j,k,l;
	double factor = 1/gamma;
	double *gridLocal = malloc(dim*sizeof(double));	
	double TermALocal, TermBLocal;
	int XCounterGlobal = 0;	

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
		gradA[i] = 0;
		gradB[i] = 0;
	}

	/* Calculate gradient for grid points */ 
	TermBLocal = 0; *TermB = 0;
	TermALocal = 0;	*TermA = 0;
	#pragma omp parallel
	{
		double *Ytmp = calloc(dim,sizeof(double));
		double stInnerMax; double stInnerCorrection = 0;
		double sum_st, sum_st_inv, tmpVal, sum_st_inv2;
		/*double *st = calloc(nH,sizeof(double)); */
		double *stInner = calloc(nH,sizeof(double));
		/*double *stInnerCheck = calloc(nH,sizeof(double));*/
		double *grad_st_private = calloc(nH*(dim+1),sizeof(double));
		double *influencePrivate = calloc(nH,sizeof(double));
		double Delta,evalTmp;
		int *idxElements = malloc(nH*sizeof(int));
		int *idxElementsBox = malloc(nH*sizeof(int));
		int numElements, numElementsBox, idxSave, idxMax = 0, sign;
		double *preCalcElems = malloc(nH*dim*MBox*sizeof(double));
		int YIdxMin, YIdxMax, idxGet, totalHyperplanes = 0;

		int XCounter=0; /* counts how much X elements were allready processed */
        double *grad_ft_private = calloc(nH*(dim+1),sizeof(double));

		/* calculate gradient for grid points */
		#pragma omp for schedule(dynamic,1) private(i,k,l) reduction(+:TermBLocal,TermALocal)
		for (j = 0; j < numBoxes; j++) {
			/* check for active hyperplanes */
			/* eval all hyperplanes for some corner point of the box */
			for (k=0; k < dim; k++) {
				Ytmp[k] = boxEvalPoints[j*3*dim + k];
			}
			stInnerMax = -DBL_MAX;
			for (i=0; i < nH; i++) {
				stInner[i] = bGamma[i] + aGamma[i]*Ytmp[0];
	
				for (k=1; k < dim; k++) {
					stInner[i] += aGamma[i+k*nH]*Ytmp[k];
				}
			}

			/* find maximum element for current grid point */
			for (i=0; i < nH; i++) {
				if (stInner[i] > stInnerMax) {
					stInnerMax = stInner[i];
					idxMax = i;
				}
			}

			if (numPointsPerBox[j+1] - numPointsPerBox[j] > 1) { 
				/* eval all hyperplanes for the point opposite to the first one */
				sign = boxEvalPoints[j*3*dim + 2*dim];
				Delta = boxEvalPoints[j*3*dim + 1*dim];
				for (i=0; i < nH; i++) {
					evalTmp = (aGamma[i]-aGamma[idxMax])*sign;
					/*stInnerCheck[i] = 0;*/
					if (evalTmp > 0) {
						stInner[i] += evalTmp*Delta;
					}
				}

				for (k=1; k < dim; k++) {
					sign = boxEvalPoints[j*3*dim + 2*dim + k];
					Delta = boxEvalPoints[j*3*dim + 1*dim + k];
					for (i=0; i < nH; i++) {
						evalTmp = (aGamma[i+k*nH]-aGamma[idxMax + k*nH])*sign;
						if (evalTmp > 0) {
							stInner[i] += evalTmp*Delta;
						}
					}
				}
			}

			/* check which hyperplanes to keep for that box */
			numElementsBox = 0;
			for (i=0; i < nH; i++) {
/*				if (stInner[i] + stInnerCheck[i] > stInnerMax + epsCalcExp) { */
				if (stInner[i]  > stInnerMax - 25) {
					idxElementsBox[numElementsBox++] = i;
				}
			}
			totalHyperplanes += numElementsBox;
			
			/* precalc elements for that box */
			for (k=0; k < dim; k++) {
				YIdxMin = INT_MAX; YIdxMax = INT_MIN;
				for (l=numPointsPerBox[j]; l < numPointsPerBox[j+1]; l++) {
					if (YIdxMin > YIdx[k + l*dim]) {
						YIdxMin = YIdx[k + l*dim];
					}
					if (YIdxMax < YIdx[k + l*dim]) {
						YIdxMax = YIdx[k + l*dim];
					}
				}

				for (l=YIdxMin; l <= YIdxMax; l++) {
					Ytmp[k] = gridLocal[k]+delta[k]*l;
					idxSave = l%MBox;
					if (k==0) {
						for (i=0; i < numElementsBox; i++) {
							preCalcElems[i + idxSave*nH] =  bGamma[idxElementsBox[i]] + aGamma[idxElementsBox[i]]*Ytmp[k]; 
						}
					} else 	{
						for (i=0; i < numElementsBox; i++) {
							preCalcElems[i + k*MBox*nH + idxSave*nH] = aGamma[idxElementsBox[i]+k*nH]*Ytmp[k];
						}
					}
				}
			}

			if (dim <= 3) {
				/* Move XCounter to the current box */
				while (XToBox[XCounter] < j) {
					XCounter++;
				}
				
				/* iterate over all samples in that box */
				while (XToBox[XCounter] == j) {
					stInnerMax = -DBL_MAX;
					for (i=0; i < numElementsBox; i++) {
						stInner[i] = bGamma[idxElementsBox[i]] + aGamma[idxElementsBox[i]]*X[XCounter];
					}
					for (k=1; k < dim; k++) {
						for (i=0; i < numElementsBox; i++) {
							stInner[i] += aGamma[idxElementsBox[i]+k*nH]*X[XCounter + (k*N)];
						}
					}

					/* find maximum element for current samples */
					for (i=0; i < numElementsBox; i++) {
						if (stInner[i] > stInnerMax) {
							stInnerMax = stInner[i];
						}
					}

					sum_st = 0; numElements = 0;
					/* calculate st only for those entries that will be non-zero */
					for (i=0; i < numElementsBox; i++) {
						if (stInner[i] - stInnerMax > -25) {
							stInner[numElements] = exp(stInner[i]-stInnerMax);
							idxElements[numElements] = i;
							sum_st += stInner[numElements++];
						}
					}

					TermALocal += XW[XCounter]*(stInnerMax + log(sum_st))*factor;
					sum_st_inv = 1/sum_st;

					/* update the gradient */
					for (i=0; i < numElements; i++) {
						idxSave = idxElementsBox[idxElements[i]];
						for (k=0; k < dim; k++) {
							grad_ft_private[idxSave + (k*nH)] += XW[XCounter]*stInner[i]*X[XCounter+(k*N)]*sum_st_inv;
						}
						grad_ft_private[idxSave + (dim*nH)] += XW[XCounter]*stInner[i]*sum_st_inv;
					}
					XCounter++;
				}
			}	

			/* iterate over all points in that box */
			for (l=numPointsPerBox[j]; l < numPointsPerBox[j+1]; l++) {
/*				printf("%d ",l); */
				stInnerMax = -DBL_MAX;
				for (k=0; k < dim; k++) {
					Ytmp[k] = gridLocal[k]+delta[k]*YIdx[k + l*dim];
				}

				idxGet = (YIdx[l*dim]%MBox)*nH;
				for (i=0; i < numElementsBox; i++) {
					//stInner[i] = bGamma[idxElementsBox[i]] + aGamma[idxElementsBox[i]]*Ytmp[0];
					stInner[i] = preCalcElems[i + idxGet];
				}
				for (k=1; k < dim; k++) {
					idxGet = (YIdx[k+l*dim]%MBox)*nH + k*MBox*nH;
					for (i=0; i < numElementsBox; i++) {
						//stInner[i] += (aGamma[idxElementsBox[i]+k*nH]*Ytmp[k]);
						stInner[i] += preCalcElems[i + idxGet];
					}
				}

				/* find maximum element for current grid point */
				for (i=0; i < numElementsBox; i++) {
					if (stInner[i] > stInnerMax) {
						stInnerMax = stInner[i];
					}
				}

				/* only calc the exponential function for elements that wont be zero afterwards */
				sum_st = 0; numElements = 0;
				for (i=0; i < numElementsBox; i++) {
					if (stInner[i] - stInnerMax > -25) {
						stInner[numElements] = exp(stInner[i]-stInnerMax);
						idxElements[numElements] = i;
						sum_st += stInner[numElements++];
					}
				}
				stInnerCorrection = exp(-stInnerMax*factor);
				tmpVal = pow(sum_st,-factor)*stInnerCorrection;
				
				TermBLocal += tmpVal; 
				sum_st_inv2 = 1/sum_st;
				sum_st_inv = tmpVal*sum_st_inv2;
				
				for (i=0; i < numElements; i++) {
					idxSave = idxElementsBox[idxElements[i]];
					influencePrivate[idxSave] += stInner[i]*sum_st_inv2;
					stInner[i] *= sum_st_inv;
					grad_st_private[idxSave] += Ytmp[0]*stInner[i];
					for (k=1; k < dim; k++) {
						grad_st_private[idxSave + (k*nH)] += Ytmp[k]*stInner[i];
					}
					grad_st_private[idxSave + dimnH] += stInner[i];
				}
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
	      	for (i=0; i < nH*(dim+1); i++) {
                gradA[i] += grad_ft_private[i];
            }

		}
		free(Ytmp); free(stInner); free(grad_st_private); free(grad_ft_private); free(influencePrivate); free(idxElements); free(idxElementsBox); free(preCalcElems); /* free(st); free(stCheck) */
	} /* end of pragma parallel */
	*TermB = TermBLocal*weight;
	for (i=0; i < nH*(dim+1); i++) {
		gradB[i] -= grad_st_tmp[i]*weight;
	}

	if (dim <= 3) {	
		/* move X pointer to elements that are not contained in any box */
		while(XToBox[XCounterGlobal] != 65535) {
			XCounterGlobal++;
		}
	}
	
	/* calculate gradient for samples X */ 
	#pragma omp parallel 
    {
   	  	double ftInnerMax;
        double sum_ft, sum_ft_inv;
        double *ft = calloc(nH,sizeof(double));
        double *ftInner = calloc(nH,sizeof(double));
        double *grad_ft_private = calloc(nH*(dim+1),sizeof(double));
		int *idxElements = malloc(nH*sizeof(int));
		int numElements, idxSave;

		/* calculate gradient for samples X */
        #pragma omp for private(i,k) schedule(dynamic) reduction(+:TermALocal)
		for (j=XCounterGlobal; j < N; j++) {
			ftInnerMax = -DBL_MAX;
			for (i=0; i < nH; i++) {
				ftInner[i] = bGamma[i] + aGamma[i]*X[j];
				
				for (k=1; k < dim; k++) {
					ftInner[i] += aGamma[i+k*nH]*X[j + (k*N)];
				}
			}

			/* find maximum element for current samples */
			for (i=0; i < nH; i++) {
				if (ftInner[i] > ftInnerMax) {
					ftInnerMax = ftInner[i];
				}
			}

			sum_ft = 0; numElements = 0;
			/* calculate ft only for those entries that will be non-zero */
			for (i=0; i < nH; i++) {
				if (ftInner[i] - ftInnerMax > -25) {
					ft[numElements] = exp(ftInner[i]-ftInnerMax);
                    idxElements[numElements] = i;
					sum_ft += ft[numElements++];
				}
			}

			TermALocal += XW[j]*(ftInnerMax + log(sum_ft))*factor;
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
	*TermA = TermALocal;
	free(grad_st_tmp); free(gridLocal); free(aGamma); free(bGamma);
}
