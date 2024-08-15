#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <headers.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#include <immintrin.h>
#include <stdint.h>

void preCondGradFloatC(int** elementList, int** elementListSize, int* numEntries, int* maxElement, int* idxEntries, float* X, float* grid, unsigned short int* YIdx, int *numPointsPerBox, float* boxEvalPoints, int numBoxes, double* a, double* b, float gamma, float weight, float* delta, int N, int M, int dim, int nH, int MBox)
{
	/* initialize elementList */
	static const int elementListIncrement = 10000000;
	free(*elementListSize); free(*elementList);
	*elementListSize = malloc(sizeof(int)); **elementListSize = elementListIncrement;
	*elementList = malloc(**elementListSize*sizeof(int));
	
	float aGamma[dim*nH], bGamma[nH];
	int i,j,k,l;
	float epsCalcExp = -250; /* this is the maximal allowed difference between the maximum and any other hyperplane in log-space; this could be possibly -700, as exp(-700) is roughly the realmin, but we are greedy and drop some hyperplanes */
    float *gridLocal = malloc(dim*sizeof(float));
    int counter = 0, savedValues = 0, savedValues2 = 0;
 
	/* initialize some variables */
    for (k=0; k < dim; k++) {
        gridLocal[k] = grid[k];
    }


	/* initialize some variables */
	for (i=0; i < nH; i++) {
		for (k=0; k < dim; k++) {
			aGamma[i+(k*nH)] = gamma*a[i+(k*nH)];
		}
		bGamma[i] = gamma*b[i];
	}

	#pragma omp parallel 
    {
   	  	float ftInnerMax;
        float *ftInner = calloc(nH,sizeof(float));
		int sizeElementList = elementListIncrement, counterLocal = 0;
		int *elementListLocal = malloc(sizeElementList*sizeof(int));
        int numElements, idxMax = 0;
		/* calculate gradient for samples */
        #pragma omp for private(i,k)
		for (j=0; j < N; j++) {
			ftInnerMax = -FLT_MAX;
			for (i=0; i < nH; i++) {
				ftInner[i] = bGamma[i] + aGamma[i]*X[j];
			}
			for (k=1; k < dim; k++) {
				for (i=0; i < nH; i++) {
					ftInner[i] += aGamma[i+k*nH]*X[j + (k*N)];
				}
			}

			/* find maximum element */
			for (i=0; i < nH; i++) {
				if (ftInner[i] > ftInnerMax) {
					ftInnerMax = ftInner[i];
					idxMax = i;
				}
			}
		
			if (sizeElementList < counterLocal + nH) {
				sizeElementList *= 2;
				elementListLocal = realloc(elementListLocal,sizeElementList*sizeof(int));
			}			

			numElements = 0;
			/* calculate ft only for those entries that are will be non-zero */
			for (i=0; i < nH; i++) {
				if (ftInner[i] - ftInnerMax > epsCalcExp) {
					elementListLocal[counterLocal++] = i;
					numElements++;
				}
			}

			maxElement[j] = idxMax;
			numEntries[j] = numElements;
		}
		
		#pragma omp critical
		{
			/* total number of elements that where added */
			counter += counterLocal;
		}

		#pragma omp barrier

		/* reallocate idxList if neccessary */
		#pragma omp single
		{
			if (**elementListSize < counter) {
				*elementList = realloc(*elementList,counter*sizeof(int));
				**elementListSize = counter;
			}
		}

#ifdef _OPENMP
        int numThreads = omp_get_num_threads();
#else
        int numThreads = 1;
#endif

		/* enforce ordered copying of memory --> keep  */
		#pragma omp for ordered schedule(static,1)
		for(j = 0; j < numThreads; j++)
		{
			#pragma omp ordered
			{
				memcpy((*elementList)+savedValues,elementListLocal,counterLocal*sizeof(int));
				savedValues += counterLocal;
			}
		}

		free(ftInner); free(elementListLocal);
	}

	#pragma omp parallel 
	{
		float stInnerMax; 
		float *stInner = calloc(nH,sizeof(float));
		int sizeElementList = elementListIncrement, counterLocal = 0;
		int *elementListLocal = malloc(sizeElementList*sizeof(int));
	   	float Delta,evalTmp;
        int *idxElements = malloc(nH*sizeof(int));
        int *idxElementsBox = malloc(nH*sizeof(int));
        int numElements, numElementsBox, idxSave, idxMax = 0, sign;
        float *preCalcElems = malloc(nH*dim*MBox*sizeof(float));
        int YIdxMin, YIdxMax, idxGet;
		float *Ytmp = calloc(dim,sizeof(float));
		int totalElements = 0;
		int *idxEntriesLocal = malloc(M*sizeof(int));
		int *maxElementLocal = malloc(M*sizeof(int));
		int *numEntriesLocal = malloc(M*sizeof(int));

		#pragma omp for schedule(dynamic,1) private(i,k,l)
        for (j=0; j < numBoxes; j++) {
            /* check for active hyperplanes */
            /* eval all hyperplanes for some corner point of the box */
            for (k=0; k < dim; k++) {
                Ytmp[k] = boxEvalPoints[j*3*dim + k];
            }
            stInnerMax = -FLT_MAX;
            for (i=0; i < nH; i++) {
                stInner[i] = bGamma[i] + aGamma[i]*Ytmp[0];
            }
            for (k=1; k < dim; k++) {
                for (i=0; i < nH; i++) {
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
					stInner[i] += evalTmp*Delta*(evalTmp > 0);
					//if (evalTmp > 0) {
					//	stInner[i] += evalTmp*Delta;
					//}
				}

				for (k=1; k < dim; k++) {
					sign = boxEvalPoints[j*3*dim + 2*dim + k];
					Delta = boxEvalPoints[j*3*dim + 1*dim + k];
					for (i=0; i < nH; i++) {
						evalTmp = (aGamma[i+k*nH]-aGamma[idxMax + k*nH])*sign;
						stInner[i] += evalTmp*Delta*(evalTmp > 0);
						//if (evalTmp > 0) {
						//	stInner[i] += evalTmp*Delta;
						//}
					}
				}
			}

            /* check which hyperplanes to keep for that box */
            numElementsBox = 0;
            for (i=0; i < nH; i++) {
                if (stInner[i] > stInnerMax + epsCalcExp) {
                    idxElementsBox[numElementsBox++] = i;
                }
            }

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
                    } else  {
                        for (i=0; i < numElementsBox; i++) {
                            preCalcElems[i + k*MBox*nH + idxSave*nH] = aGamma[idxElementsBox[i]+k*nH]*Ytmp[k];
                        }
                    }
                }
            }
       	 	
			/* iterate over all points in that box */
	        for (l=numPointsPerBox[j]; l < numPointsPerBox[j+1]; l++) {
                stInnerMax = -FLT_MAX;

                idxGet = (YIdx[l*dim]%MBox)*nH;
                for (i=0; i < numElementsBox; i++) {
                    /* stInner[i] = bGamma[idxElementsBox[i]] + aGamma[idxElementsBox[i]]*Ytmp[0];   */
                    stInner[i] = preCalcElems[i + idxGet];
                }
                for (k=1; k < dim; k++) {
                    idxGet = (YIdx[k+l*dim]%MBox)*nH + k*MBox*nH;
                    for (i=0; i < numElementsBox; i++) {
                        /* stInner[i] += (aGamma[idxElementsBox[i]+k*nH]*Ytmp[k]);  */
                        stInner[i] += preCalcElems[i + idxGet];
                    }
                }

                /* find maximum element for current grid point */
                for (i=0; i < numElementsBox; i++) {
                    if (stInner[i] > stInnerMax) {
                        stInnerMax = stInner[i];
						idxMax = i;
                    }
                }

				if (sizeElementList < counterLocal + numElementsBox) {
					sizeElementList *= 2;
					elementListLocal = realloc(elementListLocal,sizeElementList*sizeof(int));
				}

                /* only calc the exponential function for elements that wont be zero afterwards */
				numElements = 0;
				/* calculate st only for those entries that wont be zero */
				for (i=0; i < numElementsBox; i++) {
					if (stInner[i] - stInnerMax > epsCalcExp) {
						elementListLocal[counterLocal++] = idxElementsBox[i];
						numElements++;
					}
				}

				maxElementLocal[totalElements] = idxElementsBox[idxMax];
				numEntriesLocal[totalElements] = numElements;
				idxEntriesLocal[totalElements++] = l;
			}
		}
	    #pragma omp critical
        {   
            counter += counterLocal; /* final number of elements in list */
        }
  		#pragma omp barrier
 
        /* reallocate idxList if neccessary */
        #pragma omp single
        { 
         	if (**elementListSize < counter) {
                *elementList = realloc(*elementList,counter*sizeof(int));
				**elementListSize = counter;
            }
        }
#ifdef _OPENMP
        int numThreads = omp_get_num_threads();
#else
        int numThreads = 1;
#endif

        /* enforce ordered copying of memory */
        #pragma omp for ordered schedule(static,1)
        for(j = 0; j < numThreads; j++)
        {
            #pragma omp ordered
            {
				memcpy((*elementList)+savedValues,elementListLocal,counterLocal*sizeof(int));
				savedValues += counterLocal;
				memcpy(idxEntries+savedValues2,idxEntriesLocal,totalElements*sizeof(int));
				memcpy(maxElement+savedValues2+N,maxElementLocal,totalElements*sizeof(int));
				memcpy(numEntries+savedValues2+N,numEntriesLocal,totalElements*sizeof(int));
				savedValues2 += totalElements;
			}
        }

  		free(Ytmp); free(stInner); free(idxElements); free(idxElementsBox); free(preCalcElems); free(elementListLocal); free(idxEntriesLocal); free(maxElementLocal); free(numEntriesLocal);
	} /* end of pragma parallel */
	**elementListSize = counter;
}
