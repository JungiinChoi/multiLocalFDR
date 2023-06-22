#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <headers.h>
#include <limits.h>

void getMN(int dim, int sparse, int* N, int* M) {
	if (sparse) {
		if (dim==1) {
			*N = 10; *M = 5;
		} else if (dim==2) {
			*N = 10; *M = 5;
		} else if (dim==3) {
			*N = 6; *M = 4;
		} else if (dim==4) {
			*N = 4; *M = 3;
		} else if (dim==5) {
			*N = 5; *M = 2;
		} else if (dim==6) {
			*N = 4; *M = 2;
		} else if (dim==7) {
			*N = 3; *M = 2;
		} else if (dim==8) {
			*N = 2; *M = 2;
		} else if (dim==9) {
			*N = 2; *M = 2;
		} else {
			*N = 1; *M = 1;
			Rprintf("Invalid dimension for sparse grid.\n");
		}
	} else {
		if (dim==1) {
			*N = 10; *M = 50;
		} else if (dim==2) {
			*N = 10; *M = 10;
		} else if (dim==3) {
			*N = 9; *M = 5;
		} else if (dim==4) {
			*N = 6; *M = 4;
		} else if (dim==5) {
			*N = 5; *M = 3;
		} else if (dim==6) {
			*N = 4; *M = 3;
		} else if (dim==7) {
			*N = 3; *M = 3;
		} else if (dim==8) {
			*N = 3; *M = 2;
		} else if (dim==9) {
			*N = 2; *M = 2;
		} else {
			*N = 2; *M = 2;
			Rprintf("Invalid dimension for dense grid.\n"); 
		}
	}
}


void makeGridMidpoint(double* box, double* grid, int N, int dim) {
	int i,j;
	double interval;
    for (i=0; i < dim; i++) {   
        interval = (box[i+dim]-box[i])/N;
        for (j=0; j < N; j++) {
            grid[j + i*N] = box[i] + interval*(j+0.5);
        }
    }
}


void preCalcGrid(int dim, int M, int lenVec, double* vec, double* subGrid, double* preCalcGrid) {
    int i,j,k;

    for (j=0; j < dim; j++) {
        for (i=0; i < M; i++) {
            for (k=0; k < lenVec; k++) {
                preCalcGrid[i + j*M + k*dim*M] = vec[k + j*lenVec]*subGrid[i + j*M];
            }
        }
    }
}


void getSubIndex(int dim, int M, int MPow, int* subIdx) {
    int k, j, raiseIndex;

    for (j=0; j < dim; j++) {
        subIdx[j] = 0;
    }

    /* the first index always grows until M, then is reset to 0; other indices only grows if the previous one reaches M */
    for (j=1; j < MPow; j++) {
        /* increase the first one */
        raiseIndex = 0;
        if (subIdx[(j-1)*dim] < M-1) {
            subIdx[j*dim] = subIdx[(j-1)*dim]+1;
        } else {
            subIdx[j*dim] = 0;
            raiseIndex = 1;
        }
        for (k=1; k < dim; k++) {
            if (raiseIndex) {
                if (subIdx[k+(j-1)*dim] < M-1) {
                    subIdx[k+j*dim] = subIdx[k + (j-1)*dim] + 1;
                    raiseIndex = 0;
                } else {
                    subIdx[k+j*dim] = 0;
                    raiseIndex = 1;
                }
            } else {
                subIdx[k + j*dim] = subIdx[k + (j-1)*dim];
            }
        }
    }
}

void makeGridND(double *box, int N, int dim, double* sparseGrid) {
	int i, j, rep1, rep2, counter;

	double* grid = malloc(dim*(N+1)*sizeof(double));
	double delta;
	for (i=0; i < dim; i++) {
		delta = (box[i+dim] - box[i])/N;	
		// outer grid follows trapzoid rule
		for (j=0; j < N+1; j++) {
			grid[j + i*(N+1)] = box[i] + delta*j;
		}
	}

	// mimick ndgrid function from matlab
	for (i=0; i < dim; i++) {
		counter = 0;
		// repeat the arrays of the inner loop rep1 times
		for (rep1=0; rep1 < pow(N+1,dim-i-1); rep1++) {
			for (j=0; j < N+1; j++) {
				// repeat every entry in grid rep2 times
				for (rep2=0; rep2 < pow(N+1,i); rep2++) {
					sparseGrid[counter*dim + i] = grid[j + i*(N+1)];
					counter++;
				}
			}
		}
	}
	free(grid);
}

void setGridDensity(double *box, int dim, int sparseGrid, int *N, int *M, double **grid, double* weight) {
	if (sparseGrid) {
		getMN(dim,sparseGrid,N,M);
	} else {
		//if (*N==0 || *M == 0) {
		getMN(dim,sparseGrid,N,M);
		//}
	}
	*grid = malloc((*N)*(*M)*dim*sizeof(double));
	makeGridMidpoint(box, *grid, (*N)*(*M), dim); 
	*weight = 1;
	for (int i = 0; i < dim; i++) {
		*weight *= ((*grid)[N[0]*M[0]*i+1] - (*grid)[N[0]*M[0]*i]);
	}
}	

void makeGridC(double *X, unsigned short int **YIdx, unsigned short int **XToBox, int **numPointsPerBox, double **boxEvalPoints, double *ACVH, double *bCVH, double *box, int *lenY, int *numBoxes, int dim, int lenCVH, int N, int M, int NX) {

	/* ****************************************** */
	/* ********** VARIABLE DECLARATION ********** */
	/* ****************************************** */

	int i,j,k,l,numCombs,val,tmp,assign;

    int numGridPoints = pow(N+1,dim);
	int *variations = malloc(dim*numGridPoints*sizeof(int));

    int *combs = calloc(pow(2,dim)*dim,sizeof(int));
    int *idxBox = calloc(pow(2,dim),sizeof(int));
    int *factor = malloc(dim*sizeof(int));

    /* variables for pre-calculating some stuff */
    double *boxTmp = calloc(dim*2,sizeof(double));
    int *gridIdx = malloc(dim*numGridPoints*sizeof(int));
    int numPointsSubIdx;
	
	double *subGrid = malloc(dim*M*sizeof(double));
    int *subGridIdx = malloc(dim*pow(M,dim)*sizeof(int));
	
	double *sparseGrid = malloc(dim*pow(N+1,dim)*sizeof(double));
	makeGridND(box,N,dim,sparseGrid);
	double *sparseDelta = malloc(dim*sizeof(double));
	for (i=0; i < dim; i++) {
		sparseDelta[i] = (box[i+dim] - box[i])/N;
	}
    /* variables for checking boxes and adding subgrid points */
    //int maxIdx; 
	int counterNumBoxes = 0; int outsideCVH;
    int numPointsAdded = 0; int numPointsAddedOld;
    double funcEvalLocal;

	int *YBoxMin = malloc(dim*sizeof(int));
	int *YBoxMax = malloc(dim*sizeof(int));
	double *boxMax = malloc(dim*sizeof(double));
	double *boxMaxOuter = malloc(dim*sizeof(double));
	double *boxMinOuter = malloc(dim*sizeof(double));
	double *boxMin = malloc(dim*sizeof(double));
	double epsilon = pow(10,-10);
	double yTmp[dim];
    /* create arrays that are supposed to be returned to the calling function */
	*numBoxes = pow(N,dim);
	*lenY = *numBoxes*pow(M,dim);
	*YIdx = malloc(*lenY*dim*sizeof(unsigned short int));
    *XToBox = malloc(NX*sizeof(unsigned short int));
    unsigned short int *XToBoxOuter = malloc(NX*sizeof(unsigned short int));
	for (i=0; i < NX; i++) {
		XToBoxOuter[i] = USHRT_MAX;
		(*XToBox)[i] = USHRT_MAX;
	}
	
#ifdef __AVX__
	int nB = ((int) (*numBoxes/8) + 1)*8;
#else
	int nB = *numBoxes;
#endif
    *numPointsPerBox = calloc(*numBoxes+1,sizeof(int));
	*boxEvalPoints = calloc(*numBoxes*3*dim,sizeof(double)); 

    /* ********************************************************************** */
    /* *********** PRECALCULATE STUFF --> SPEEDS UP MAIN FOR LOOP *********** */
    /* ********************************************************************** */

    /* create all combinations of a binary vector of length dim: 000, 001, 010, 011, ... */
    numCombs = pow(2,dim);
    for (i=1; i < numCombs; i++) {
        val = i;
        for (j=dim-1; j >= 0; j--) {
            tmp = pow(2,j);
            if (tmp <= val) {
                combs[i*dim+j] = 1;
                val = val-tmp;
            }
        }
    }

    for (j=0; j < dim; j++) {
        factor[j] = pow(N+1,j);
    }

    for (i=0; i < numCombs; i++) {
        for(j=0; j < dim; j++) {
            idxBox[i] += combs[i*dim+j]*factor[j];
        }
    }
    free(factor); free(combs);

    /* make subgrid for midpoint rule (grid inside each box) */
    for (j=0; j < dim; j++) {
        boxTmp[dim+j] = (box[j+dim] - box[j])/N;
   	}
    makeGridMidpoint(boxTmp,subGrid,M,dim);
    free(boxTmp);

    /* index for selecting the subgrid */
    numPointsSubIdx = pow(M,dim);
    getSubIndex(dim,M,numPointsSubIdx,subGridIdx);
    getSubIndex(dim,(N+1),numGridPoints,gridIdx);

    /* ************************************************************************************ */
    /* ************* MAIN LOOP: ADDS ALL POINTS LYING INSIDE THE CONVEX HULL ************** */
    /* ************************************************************************************ */

	counterNumBoxes = 0; numPointsAdded = 0;
   	for (l=0; l < numGridPoints; l++) {
        /* boxes can't correspond to grid points that lie on the right boundary */
      /*  for (i=0; i < dim; i++) {
            if (variations[i+dim*l]>maxIdx) {
                maxIdx = variations[i+dim*l];
            }
        }
        if (maxIdx==N) {
			printf("test");
            continue;
        }*/

		/* reset box min and max */
		for (j=0; j < dim; j++) {
			YBoxMax[j] = 0; YBoxMin[j] = M-1;
		}

        numPointsAddedOld = numPointsAdded;
        /* check points in subgrid */
        for (j=0; j < numPointsSubIdx; j++) {
            outsideCVH = 0;
			for (k=0; k < dim; k++) {
				yTmp[k] = sparseGrid[l*dim + k] + subGrid[subGridIdx[k + j*dim] + k*M];
			}
			/* check each point */
            for (i=0; i < lenCVH; i++) {
                funcEvalLocal = -bCVH[i];
                for (k=0; k < dim; k++) {
                    funcEvalLocal += ACVH[i+k*lenCVH]*yTmp[k];
                }
                if (funcEvalLocal > 0) {
                    outsideCVH = 1;
					/*if (l==0 && j == 0) {
						Rprintf("Failed for %d with %.4e\n",i,funcEvalLocal);
					}*/
                    break; /* as soon as one violating hyperplane is found, we can stop the search and continue with the next grid point */
                }
            }
            /* if inside convex hull, add point to list of points */
            if (outsideCVH==0) {
                for (k=0; k < dim; k++) {
/*                  (*Y)[k + numPointsAdded*dim] = sparseGrid[l*dim + k] + subGrid[subGridIdx[k + j*dim] + k*N]; */
                    (*YIdx)[k + numPointsAdded*dim] = subGridIdx[k + j*dim] + gridIdx[k + l*dim]*M;
					/* check for min and max */
					if (YBoxMin[k] > subGridIdx[k+j*dim]) {
						YBoxMin[k] = subGridIdx[k+j*dim];
					}
					if (YBoxMax[k] < subGridIdx[k+j*dim]) {
						YBoxMax[k] = subGridIdx[k+j*dim];
					}
                }
                /* save subGridIdx (1-D index from 1 to M^d) */
                numPointsAdded++;
				//Rprintf("%d: %d, %d added\n",numPointsAdded, l, j);
		   	}
		}
		/* find box max min of outer box boundaries */
		for (j=0; j < dim; j++) {
			boxMinOuter[j] = sparseGrid[l*dim+j];
			boxMaxOuter[j] = sparseGrid[l*dim+j] + sparseDelta[j];
		}

		/* find box max and min, based on the sparse grid defining box boundaries (and not on the subgrid inside the box) 
		 * With that definition, we can use the estimates of active hyperplanes for each box also when evaluating X */
		for (j=0; j < dim; j++) {
			boxMin[j] = sparseGrid[l*dim + j]+subGrid[YBoxMin[j] + j*M]; 
			boxMax[j] = sparseGrid[l*dim + j]+subGrid[YBoxMax[j] + j*M];
		}

        /* only add the box, if at least one grid point was added */
        if (numPointsAdded > numPointsAddedOld) {
			/* iterate over all X values and check whether they belong to the current box */
			for (j=0; j < NX; j++) {
				/* if still unassigned */
				if (XToBoxOuter[j] == USHRT_MAX) {
					assign = 1;
					for (k=0; k < dim; k++) {
						if (X[k*NX + j] > boxMaxOuter[k]+epsilon || X[k*NX+j] < boxMinOuter[k]-epsilon) {
							assign = 0;
							break;
						}
					}
					if (assign==1) {
						XToBoxOuter[j] = counterNumBoxes;
						assign = 1;
						/* check for inner box */
						for (k=0; k < dim; k++) {
							if (X[k*NX + j] > boxMax[k]+epsilon || X[k*NX+j] < boxMin[k]-epsilon) {
								assign = 0;
								break;
							}
						}
						if (assign==1) {
							(*XToBox)[j] = counterNumBoxes;
						}
					}
				}
			}
	    	(*numPointsPerBox)[counterNumBoxes+1] = numPointsAdded;
			for (j=0; j < dim; j++) {
				/* choose the point closest to zero, as it most probably lies inside the convex hull of X */
				/* for each box three points: The point to evaluate, the delta between the two corners and the sign indicating in which corner we start */
				if (fabs(boxMin[j]) < fabs(boxMax[j])) {
					(*boxEvalPoints)[counterNumBoxes*3*dim + j] = boxMin[j];
					(*boxEvalPoints)[counterNumBoxes*3*dim + 2*dim + j] = 1;
				} else {
					(*boxEvalPoints)[counterNumBoxes*3*dim + j] = boxMax[j];
					(*boxEvalPoints)[counterNumBoxes*3*dim + 2*dim + j] = -1;
				}
			    (*boxEvalPoints)[counterNumBoxes*3*dim + dim + j] = boxMax[j]-boxMin[j];
			}
			counterNumBoxes++;
 		}
    }
	// resize arrays to the final number of grid points/boxes
    if (counterNumBoxes != *numBoxes) {
#ifdef __AVX__
		nB = ((int) (counterNumBoxes/8) + 1)*8;
#else
		nB = counterNumBoxes;
#endif
        *numPointsPerBox = realloc(*numPointsPerBox,(nB+1)*sizeof(int));
		*boxEvalPoints = realloc(*boxEvalPoints,nB*dim*3*sizeof(double)); 
    }
    *numBoxes = counterNumBoxes;

    if (numPointsAdded != *lenY) {
        *YIdx = realloc(*YIdx,numPointsAdded*dim*sizeof(unsigned short int));
    }
    *lenY = numPointsAdded;

    free(idxBox);
    free(variations);
    free(gridIdx);
	free(YBoxMax); free(YBoxMin);
	free(sparseGrid);
	free(subGrid); free(subGridIdx);
	free(boxMax); free(boxMin); free(boxMaxOuter); free(boxMinOuter);
	free(XToBoxOuter);
	free(sparseDelta);
}
