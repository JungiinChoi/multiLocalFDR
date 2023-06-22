#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <limits.h>
#include <sys/time.h>
#include <string.h>
#include <stdint.h>
#include <R.h>
#include <headers.h>

#ifdef __AVX__
#include "avx_mathfun.h"
#define ALIGN 32

#define _mm256_full_hadd_ps(v0, v1) \
        _mm256_hadd_ps(_mm256_permute2f128_ps(v0, v1, 0x20), \
                       _mm256_permute2f128_ps(v0, v1, 0x31))

void evalHyperplaneScalar(float* aGamma, float* bGamma, float* ftInner, float* X, float* XW, float* grad_ft_private, float* TermALocal, int* idxElements, int j, int nHLocal, int nH, int* idxElementsBox, int dim, int N, float factor) {
	int i,k,idxSave,numElements;
	float sum_ft_scalar, sum_ft_scalar_inv, tmpVal, ftInnerMax;
	ftInnerMax = -FLT_MAX;
	numElements = 0;
	for (i=0; i < nHLocal; i++) {
		tmpVal = bGamma[i];
		for (k=0; k < dim; k++) {
			tmpVal += aGamma[i*dim+k]*X[j + (k*N)];
		}
		if (tmpVal - (ftInnerMax) > -25) {
			if (tmpVal > ftInnerMax) {
				ftInnerMax = tmpVal;
			}
			ftInner[numElements] = tmpVal;
			idxElements[numElements++] = i;
		}
	}

	sum_ft_scalar = 0;	
	// calculate ft only for those entries that will be non-zero
	for (i=0; i < numElements; i++) {
		ftInner[i] = expf(ftInner[i]-ftInnerMax);
		sum_ft_scalar += ftInner[i];
	}

	*TermALocal += XW[j]*(ftInnerMax + logf(sum_ft_scalar))*factor;
	sum_ft_scalar_inv = 1/sum_ft_scalar*XW[j];

	// update the gradient
	if (nH!=nHLocal) {
		for (i=0; i < numElements; i++) {
			idxSave = idxElementsBox[idxElements[i]];
			for (k=0; k < dim; k++) {
				grad_ft_private[idxSave + (k*nH)] += ftInner[i]*X[j+(k*N)]*sum_ft_scalar_inv;
			}
			grad_ft_private[idxSave + (dim*nH)] += ftInner[i]*sum_ft_scalar_inv;
		}
	} else	{
		for (i=0; i < numElements; i++) {
            idxSave = idxElements[i];
            for (k=0; k < dim; k++) {
                grad_ft_private[idxSave + (k*nH)] += ftInner[i]*X[j+(k*N)]*sum_ft_scalar_inv;
            }
            grad_ft_private[idxSave + (dim*nH)] += ftInner[i]*sum_ft_scalar_inv;
        }
	}
}

void evalHyperplane(float* aGamma, float* bGamma, int* numElements, int* idxElements, float* ftInner, float* XAligned, int dim, int nH, int N, __m256* sum_ft, __m256* ftMax) {
    int idxA, i, k;
    __m256 ft, cmp, a, x[dim], val1;
    __m256 delta = _mm256_set1_ps(-25);

    *numElements = 0;
    *ftMax = _mm256_set1_ps(-FLT_MAX);
	for (k=0; k < dim; k++) {
		x[k] = _mm256_loadu_ps(XAligned + N*k);
	}
    for (i=0; i < nH; i++) {
        idxA = i*dim;
        ft = _mm256_set1_ps(*(bGamma + i));
        // ftTmp = bGamma[i] + aGamma[idxA]*X[idxB];
        for (k=0; k < dim; k++) {
            //x = _mm256_loadu_ps(XAligned + N*k); // pull values from X for one dimension and 8 values
            a = _mm256_set1_ps(*(aGamma + idxA + k)); // fill a with the same scalar
#ifdef __AVX2__
            ft = _mm256_fmadd_ps(x[k],a,ft); // combined multiply+add
#else
			ft = _mm256_add_ps(ft,_mm256_mul_ps(x[k],a));
#endif
        }

        // ftTmp > ftInnerMax - 25
        val1 = _mm256_add_ps(*ftMax,delta);
        cmp =_mm256_cmp_ps(ft,val1,_CMP_GT_OQ);

        // check if any value surpasses the maximum
        if (_mm256_movemask_ps(cmp) > 0) {
            // update max vals
            *ftMax = _mm256_max_ps(ft,*ftMax);
            // save ft values for later use
            _mm256_store_ps(ftInner + 8*(*numElements), ft);
            idxElements[(*numElements)++] = i;
	   	}
    }
    // set all values to zero
    *sum_ft = _mm256_setzero_ps();
    // calculate exp(ft) and sum_ft
    for (i=0; i < *numElements; i++) {
        ft = _mm256_load_ps(ftInner + 8*i);
        ft = exp256_ps(_mm256_sub_ps(ft,*ftMax)); // ftInner[i] = exp(ftInner[i]-ftInnerMax);
        _mm256_store_ps(ftInner + 8*i,ft);
        *sum_ft = _mm256_add_ps(*sum_ft,ft); // sum_ft += ftInner[i];
    }
}

void evalHyperplaneY(float* aGamma, float* bGamma, int* numElements, int* idxElements, float* ftInner, float* XAligned, int dim, int nH, int N, __m256* sum_ft, __m256* ftMax, int* elemsIdx) {
    int idxA, i, k;
    __m256 ft, cmp, a, x[dim], val1;
    __m256 delta = _mm256_set1_ps(-25);

    *numElements = 0;
    *ftMax = _mm256_set1_ps(-FLT_MAX);
	for (k=0; k < dim; k++) {
		x[k] = _mm256_loadu_ps(XAligned + N*k);
	}
    for (i=0; i < nH; i++) {
        idxA = i*dim;
        ft = _mm256_set1_ps(*(bGamma + i));
        // ftTmp = bGamma[i] + aGamma[idxA]*X[idxB];
        for (k=0; k < dim; k++) {
            a = _mm256_set1_ps(*(aGamma + idxA + k)); // fill a with the same scalar
#ifdef __AVX2__
            ft = _mm256_fmadd_ps(x[k],a,ft); // combined multiply+add
#else
			ft = _mm256_add_ps(ft,_mm256_mul_ps(x[k],a));
#endif
        }

        // ftTmp > ftInnerMax - 25
        val1 = _mm256_add_ps(*ftMax,delta);
        cmp =_mm256_cmp_ps(ft,val1,_CMP_GT_OQ);

        // check if any value surpasses the maximum
        if (_mm256_movemask_ps(cmp) > 0) {
            // update max vals
            *ftMax = _mm256_max_ps(ft,*ftMax);
            // save ft values for later use
            _mm256_store_ps(ftInner + 8*(*numElements), ft);
            idxElements[(*numElements)++] = elemsIdx[i];
        }
    }
    // set all values to zero
    *sum_ft = _mm256_setzero_ps();
    // calculate exp(ft) and sum_ft
    for (i=0; i < *numElements; i++) {
        ft = _mm256_load_ps(ftInner + 8*i);
        ft = exp256_ps(_mm256_sub_ps(ft,*ftMax)); // ftInner[i] = exp(ftInner[i]-ftInnerMax);
        _mm256_store_ps(ftInner + 8*i,ft);
        *sum_ft = _mm256_add_ps(*sum_ft,ft); // sum_ft += ftInner[i];
    }
}

void calcGradient(int numElements, int* idxElements, float* ftInner, float* XAligned, float* grad_ft_private, int dim, int nH, int N, __m256 sum_ft_inv) {
    float *t;
    __m256 val1,ft,x[dim];
    int i,k,idxSave;
	for (k=0; k < dim; k++) {
		x[k] = _mm256_loadu_ps(XAligned + N*k);
	}
    for (i=0; i < numElements; i++) {
        idxSave = idxElements[i];
        ft = _mm256_load_ps(ftInner + 8*i);
        ft = _mm256_mul_ps(ft,sum_ft_inv);
        for (k=0; k < dim; k++) {
            val1 = _mm256_mul_ps(ft,x[k]);
            t = (float*) &val1;
            grad_ft_private[idxSave + (k*nH)] += t[0] + t[1] + t[2] + t[3] + t[4] + t[5] + t[6] + t[7];
        }
        t = (float*) &ft;
        grad_ft_private[idxSave + (dim*nH)] += t[0] + t[1] + t[2] + t[3] + t[4] + t[5] + t[6] + t[7];
    }
}


void calcInfluence(int numElements, int* idxElements, float* ftInner, float* influencePrivate, __m256 sum_ft_inv2) {
    float* t;
    __m256 ft;
    int i, idxSave;
    for (i=0; i < numElements; i++) {
        idxSave = idxElements[i];
        ft = _mm256_load_ps(ftInner + 8*i);
        ft = _mm256_mul_ps(ft,sum_ft_inv2); //st[i]*sum_st_inv2
        t = (float*) &ft;
        influencePrivate[idxSave] += t[0] + t[1] + t[2] + t[3] + t[4] + t[5] + t[6] + t[7];
    }
}

void findMaxVal(float* aGamma, float* bGamma, float* ftInner, float* X, int dim, int nH, int N, __m256* ftMax, float* boxEvalPoints, int* numPointsPerBox, int* idxElementsBox, int* numElementsBox) {
  	int idxA, i, k;
    __m256 ft, cmp, a, x, val1;

	float idxMax[8];
	int mask;
    *ftMax = _mm256_set1_ps(-FLT_MAX);
	memset(idxMax,0,8*sizeof(float));
	memset(numElementsBox,0,8*sizeof(int));
	for (i=0; i < nH; i++) {
        idxA = i*dim;
        ft = _mm256_set1_ps(*(bGamma + i));
        // ftTmp = bGamma[i] + aGamma[idxA]*X[idxB];
        for (k=0; k < dim; k++) {
            x = _mm256_loadu_ps(X + N*k); // pull values from X for one dimension and 8 values
            a = _mm256_set1_ps(*(aGamma + idxA + k)); // fill a with the same scalar
#ifdef __AVX2__
            ft = _mm256_fmadd_ps(x,a,ft); // combined multiply+add
#else
			ft = _mm256_add_ps(ft,_mm256_mul_ps(x,a));
#endif
        }
		_mm256_store_ps(ftInner+8*i,ft);
		cmp =_mm256_cmp_ps(ft,*ftMax,_CMP_GT_OQ);
//		_mm256_maskstore_epi32(idxMax,_mm256_castps_si256(cmp),_mm256_set1_epi32(i)); AVX2 only
		_mm256_maskstore_ps(idxMax,_mm256_castps_si256(cmp),_mm256_set1_ps(i));
		*ftMax = _mm256_max_ps(ft,*ftMax);
	}
	int idxMaxCast[8];
	for (i=0; i < 8; i++) {
		idxMaxCast[i] = (int) idxMax[i];
	}

	__m256 sign[dim],numPoints,Delta[dim],aGammaMax[dim],evalTmp,zeros,ones;
	for (k=0; k < dim; k++) {
		aGammaMax[k] = _mm256_set_ps(aGamma[idxMaxCast[7]*dim+k],aGamma[idxMaxCast[6]*dim+k],aGamma[idxMaxCast[5]*dim+k],aGamma[idxMaxCast[4]*dim+k],aGamma[idxMaxCast[3]*dim+k],aGamma[idxMaxCast[2]*dim+k],aGamma[idxMaxCast[1]*dim+k],aGamma[idxMaxCast[0]*dim+k]);
		sign[k] = _mm256_set_ps(boxEvalPoints[7*3*dim + 2*dim + k],boxEvalPoints[6*3*dim + 2*dim + k],boxEvalPoints[5*3*dim + 2*dim + k],boxEvalPoints[4*3*dim + 2*dim + k],boxEvalPoints[3*3*dim + 2*dim + k],boxEvalPoints[2*3*dim + 2*dim + k],boxEvalPoints[1*3*dim + 2*dim + k],boxEvalPoints[0*3*dim + 2*dim + k]);
		Delta[k] = _mm256_set_ps(boxEvalPoints[7*3*dim + 1*dim + k],boxEvalPoints[6*3*dim + 1*dim + k],boxEvalPoints[5*3*dim + 1*dim + k],boxEvalPoints[4*3*dim + 1*dim + k],boxEvalPoints[3*3*dim + 1*dim + k],boxEvalPoints[2*3*dim + 1*dim + k],boxEvalPoints[1*3*dim + 1*dim + k],boxEvalPoints[0*3*dim + 1*dim + k]);
	}
	zeros = _mm256_setzero_ps();
	ones = _mm256_set1_ps(1);
	numPoints = _mm256_set_ps(numPointsPerBox[8]-numPointsPerBox[7]>1,numPointsPerBox[7]-numPointsPerBox[6]>1,numPointsPerBox[6]-numPointsPerBox[5]>1,numPointsPerBox[5]-numPointsPerBox[4]>1,numPointsPerBox[4]-numPointsPerBox[3]>1,numPointsPerBox[3]-numPointsPerBox[2]>1,numPointsPerBox[2]-numPointsPerBox[1]>1,numPointsPerBox[1]-numPointsPerBox[0]>1);

	*ftMax = _mm256_add_ps(*ftMax,_mm256_set1_ps(-25));
	for (i=0; i < nH; i++) {
		ft = _mm256_load_ps(ftInner+8*i);
		for (k=0; k < dim; k++) {
			a = _mm256_set1_ps(*(aGamma + i*dim + k));
			evalTmp = _mm256_mul_ps(_mm256_sub_ps(a,aGammaMax[k]),sign[k]);
			cmp = _mm256_blendv_ps(zeros,ones,_mm256_cmp_ps(evalTmp,zeros,_CMP_GT_OQ));
			val1 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(evalTmp,Delta[k]),cmp),numPoints);
			ft = _mm256_add_ps(ft,val1); 
		}
//		_mm256_store_ps(ftInner+8*i,ft);
	    cmp =_mm256_cmp_ps(ft,*ftMax,_CMP_GT_OQ);
		mask = _mm256_movemask_ps(cmp);
		if (mask !=0) {
			float* t = (float*) &cmp;
			for (int l = 0; l < 8; l++) {
				if (((int) t[l]) !=0) {
					idxElementsBox[l*nH + numElementsBox[l]++] = i;
				}
			}
		}
	}
}

void calcGradAVXC(double* gradA, double* gradB, double* influence, double* TermA, double* TermB, float* X, float* XW, float* grid, unsigned short int* YIdx, int *numPointsPerBox, float* boxEvalPoints, unsigned short int *XToBox, int numBoxes, double* a, double* b, float gamma, float weight, float* delta, int N, int M, int dim, int nH)
{
	float *grad_st_tmp = calloc(nH*(dim+1),sizeof(float));
	float *aGamma = malloc(dim*nH*sizeof(float)); 
	float *bGamma = malloc(nH*sizeof(float));
	int dimnH = dim*nH;
	int i,j,k;
	float factor = 1/gamma;
	float *gridLocal = malloc(dim*sizeof(float));	
	float TermALocal, TermBLocal;
	int XCounterGlobal = 0;	

	// initialize some variables
	for (k=0; k < dim; k++) {
		gridLocal[k] = grid[k];
	}

	for (i=0; i < nH; i++) {
		for (j=0; j < dim; j++) {
			aGamma[i*dim + j] = (float) gamma*a[i*dim+j];
		}
		bGamma[i] = (float) gamma*b[i];
		influence[i] = 0;
	}

	// Calculate gradient for grid points
	TermBLocal = 0; *TermB = 0;
	TermALocal = 0;	*TermA = 0;
	long int totalHyperplanes;
	int countA, countB, countC, countD; countA = countB = countC = countD = 0;
	#pragma omp parallel
	{
		float *Ytmp = malloc(8*dim*sizeof(float));
		float stInnerMax; float stInnerCorrection = 0;
		float sum_st, sum_st_inv, tmpVal, sum_st_inv2;
		float *stInner;
		float *grad_st_private = calloc(nH*(dim+1),sizeof(float));
		float *influencePrivate = calloc(nH,sizeof(float));
		//float Delta;
        float *aLocal, *bLocal;
		alloc_aligned_mem(8*nH,ALIGN,&stInner); 
		alloc_aligned_mem(dim*nH,ALIGN,&aLocal); 
		alloc_aligned_mem(8*nH,ALIGN,&bLocal); 
		int *idxElements = malloc(8*nH*sizeof(int));
		int *idxElementsBox = malloc(8*nH*sizeof(int));
		int numElements, numElementsBox[8], idxSave; //idxMax;
		//int YIdxMin, YIdxMax, idxGet;
		int s1, s2, numPoints;
		int l,n,m;
		// variables AVX
		float *t;
    	__m256 stMax_,sum_st_,xw,factorNeg_,val1,val2,sum_st_inv_,sum_st_inv2_,ones,stInnerCorrection_,factor_;
        factorNeg_ = _mm256_set1_ps(-factor);
     	factor_ = _mm256_set1_ps(factor);
        ones = _mm256_set1_ps(1);

		int XCounter=0; // counts how much X elements were allready processed
        float *grad_ft_private = calloc(nH*(dim+1),sizeof(float));

		// calculate gradient for grid points
		#pragma omp for schedule(dynamic,1) private(j,i,k) reduction(+:TermBLocal,TermALocal,totalHyperplanes,countA,countB,countC,countD)
		for (n = 0; n < numBoxes; n+=8) {
			// check for active hyperplanes 
			// eval all hyperplanes for some corner point of the box
			memset(Ytmp, 0, 8*dim*sizeof(float));
			for (m = 0; m < fmin(n+8,numBoxes)-n; m++) {
				for (k=0; k < dim; k++) {
					Ytmp[m + k*8] = boxEvalPoints[(n+m)*3*dim + k];
				}
			}
			findMaxVal(aGamma, bGamma, stInner, Ytmp, dim, nH, 8, &stMax_,boxEvalPoints+n*3*dim,numPointsPerBox+n,idxElementsBox,numElementsBox);

			for (m = 0; m < fmin(n+8,numBoxes)-n; m++) {
				j = n+m;
				
				totalHyperplanes += numElementsBox[m]*(numPointsPerBox[j+1] - numPointsPerBox[j]);
				//storeHyperplanes(aTransLocal,bLocal,bGamma,aTransGamma,idxElementsBox,numElementsBox,m,dim,nH);
				// save hyperplanes in one compact vector
				for (i=0; i < numElementsBox[m]; i++) {
					bLocal[i] = bGamma[idxElementsBox[i+m*nH]];
					for (k=0; k < dim; k++) {
						aLocal[i*dim + k] = aGamma[idxElementsBox[i+m*nH]*dim + k];
					}
				}

				// Move XCounter to the current box
				if (dim <= 3) {
					while (XToBox[XCounter] < j) {
						XCounter++;
					}
			
					// iterate over all samples in that box 
					while (N > XCounter+7 && XToBox[XCounter+7]==j) {
						countA++;
						evalHyperplaneY(aLocal, bLocal, &numElements, idxElements, stInner, X+XCounter, dim, numElementsBox[m], N, &sum_st_, &stMax_,idxElementsBox+m*nH);
						xw = _mm256_loadu_ps(XW + XCounter);
						val1 = _mm256_add_ps(stMax_,log256_ps(sum_st_)); //(ftInnerMax + log(sum_ft))
						val2 = _mm256_mul_ps(_mm256_mul_ps(xw,val1),factor_); //*XW[j]*factor
						t = (float*) &val2;
						TermALocal += t[0] + t[1] + t[2] + t[3] + t[4] + t[5] + t[6] + t[7];

						sum_st_inv_ = _mm256_div_ps(xw,sum_st_); // sum_st_inv = 1/sum_st*XW[j]

						calcGradient(numElements, idxElements, stInner, X + XCounter, grad_ft_private, dim, nH, N, sum_st_inv_);
						XCounter+=8;
					}
				
					while (XToBox[XCounter] == j) {
						countB++;
						evalHyperplaneScalar(aLocal, bLocal, stInner, X, XW, grad_ft_private, &TermALocal, idxElements, XCounter, numElementsBox[m], nH, idxElementsBox+m*nH, dim, N, factor);
						XCounter++;
					}
				}

				s1 = numPointsPerBox[j]; s2 = numPointsPerBox[j+1]; numPoints = s2 - s1;
				// iterate over all points in that box
				for (l = 0; l < numPoints - numPoints%8; l+=8) {
					countC++;
					for (k=0; k < dim; k++) {
						for (i=0; i < 8; i++) {
						   Ytmp[i + k*8] = gridLocal[k]+delta[k]*YIdx[k + (s1 + l + i)*dim];
						}
					}
					evalHyperplaneY(aLocal, bLocal, &numElements, idxElements, stInner, Ytmp, dim, numElementsBox[m], 8, &sum_st_, &stMax_,idxElementsBox+m*nH);
					stInnerCorrection_ = exp256_ps(_mm256_mul_ps(stMax_,factorNeg_)); // stInnerCorrection = exp(stInnerMax*-factor);
					val1 = _mm256_mul_ps(exp256_ps(_mm256_mul_ps(log256_ps(sum_st_),factorNeg_)),stInnerCorrection_);  // tmpVal = pow(sum_st,-factor)*stInnerCorrection;

					t = (float *) &val1;
					TermBLocal += t[0] + t[1] + t[2] + t[3] + t[4] + t[5] + t[6] + t[7]; 
					//_mm256_storeu_ps(evalFunc+l+s1,val1); // evalGrid[j] = tmpVal;
					sum_st_inv2_ = _mm256_div_ps(ones,sum_st_); // sum_st_inv = 1/sum_st;
					sum_st_inv_ = _mm256_mul_ps(val1,sum_st_inv2_); // sum_st_inv2 = tmpVal*sum_st_inv;

					calcGradient(numElements, idxElements, stInner, Ytmp, grad_st_private, dim, nH, 8, sum_st_inv_);
					calcInfluence(numElements, idxElements, stInner, influencePrivate, sum_st_inv2_);
				}

				for (l = l; l < numPoints; l++) { 
					countD++;
					stInnerMax = -FLT_MAX;
					for (k=0; k < dim; k++) {
						Ytmp[k] = gridLocal[k]+delta[k]*YIdx[k + (l+s1)*dim];
					}

					numElements = 0;
					for (i=0; i < numElementsBox[m]; i++) {
						tmpVal = bLocal[i];
					   	for (k=0; k < dim; k++) {
							tmpVal += aLocal[i*dim + k]*Ytmp[k];
						}
						if (tmpVal - stInnerMax > -25) {
							if (tmpVal > stInnerMax) {
								stInnerMax = tmpVal;
							}
							stInner[numElements] = tmpVal;
							idxElements[numElements++] = i;
						}
					}				

					// only calc the exponential function for elements that wont be zero afterwards
					sum_st = 0;
					for (i=0; i < numElements; i++) {
						stInner[i] = expf(stInner[i]-stInnerMax);
						sum_st += stInner[i];
					}
					stInnerCorrection = expf(-stInnerMax*factor);
					tmpVal = pow(sum_st,-factor)*stInnerCorrection;
					
					TermBLocal += tmpVal; //evalFunc[l+s1] = tmpVal;
					sum_st_inv2 = 1/sum_st;
					sum_st_inv = tmpVal*sum_st_inv2;
					
					for (i=0; i < numElements; i++) {
						idxSave = idxElementsBox[idxElements[i]+m*nH];
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
		}
		#pragma omp critical
		{
			for (i=0; i < nH; i++) {
				influence[i] += (double) influencePrivate[i];
			}
			for (i=0; i < nH*(dim+1); i++) {
				grad_st_tmp[i] += grad_st_private[i];
			}
			for (i=0; i < nH*(dim+1); i++) {
				gradA[i] += (double) grad_ft_private[i];
			}

		}
		free(Ytmp); free(grad_st_private); free(grad_ft_private); free(influencePrivate); free(idxElements); free(idxElementsBox); 
		free_aligned_mem(aLocal); free_aligned_mem(bLocal); free_aligned_mem(stInner);  
	} // end of pragma parallel 
	*TermB = (double) TermBLocal*weight;
	for (i=0; i < nH*(dim+1); i++) {
		gradB[i] -= (double) grad_st_tmp[i]*weight;
	}

	// move X pointer to elements that are not contained in any box 
	if (dim <= 3) {
		while(XToBox[XCounterGlobal] != 65535) {
			XCounterGlobal++;
		}
	} else {
		XCounterGlobal = 0;
	}
	int XCounterGlobalStart = XCounterGlobal;
	
	// calculate gradient for samples X 
	#pragma omp parallel
    {
        float *grad_ft_private = calloc(nH*(dim+1),sizeof(float));
		float *ftInner;
		alloc_aligned_mem(8*nH,ALIGN,&ftInner); 
		int *idxElements = malloc(nH*sizeof(int));
        int numElements;
        float *t;
		int *idxElementsBox = NULL;
        __m256 ftMax,sum_ft,xw,factor_,val1,val2,sum_ft_inv;
        factor_ = _mm256_set1_ps(factor);
  		#pragma omp for schedule(dynamic) reduction(+:TermALocal)
        for (j=XCounterGlobalStart; j < N-8; j+=8) {
			#pragma omp atomic
			XCounterGlobal+=8;
            evalHyperplane(aGamma, bGamma, &numElements, idxElements, ftInner, X + j, dim, nH, N, &sum_ft, &ftMax);
            xw = _mm256_loadu_ps(XW + j);
            val1 = _mm256_add_ps(ftMax,log256_ps(sum_ft)); //(ftInnerMax + log(sum_ft))
            val2 = _mm256_mul_ps(_mm256_mul_ps(xw,val1),factor_); // *XW[j]*factor
            t = (float*) &val2;
            TermALocal += t[0] + t[1] + t[2] + t[3] + t[4] + t[5] + t[6] + t[7];
            sum_ft_inv = _mm256_div_ps(xw,sum_ft); // sum_ft_inv = 1/sum_ft*XW[j]

            calcGradient(numElements, idxElements, ftInner, X + j, grad_ft_private, dim, nH, N, sum_ft_inv);
        }
		#pragma omp barrier
        #pragma omp for reduction(+:TermALocal)
        for (j=XCounterGlobal; j < N; j++) {
			evalHyperplaneScalar(aGamma, bGamma, ftInner, X, XW, grad_ft_private, &TermALocal, idxElements, j, nH, nH, idxElementsBox, dim, N, factor);
		}
		#pragma omp critical
        {
			for (i=0; i < nH*(dim+1); i++) {
                gradA[i] += (double) grad_ft_private[i];
            }
        }
		free_aligned_mem(ftInner); free(grad_ft_private); free(idxElements);
	}
	*TermA = (double) TermALocal;
	free(grad_st_tmp); free(gridLocal); free(aGamma); free(bGamma);
}
#else

void calcGradAVXC(double* gradA, double* gradB, double* influence, double* TermA, double* TermB, float* X, float* XW, float* grid, unsigned short int* YIdx, int *numPointsPerBox, float* boxEvalPoints, unsigned short int *XToBox, int numBoxes, double* a, double* b, float gamma, float weight, float* delta, int N, int M, int dim, int nH){
	// empty function --> AVX not used
}
#endif
