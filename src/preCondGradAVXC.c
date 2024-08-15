#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#include <stdint.h>
#include <sys/time.h>

#ifdef __AVX__
#include <headers.h>

#define ALIGN 32

struct maxS  {
    float val;
    unsigned short int idx;
};

// horizontal max function; returns index and value in the struct maxS
struct maxS avx_max_struct(__m256 a1,__m256 ft) {
    __m256 a2,a3,a4,b1,b2,c1;
    float* t;
    struct maxS p;
    a2 = _mm256_shuffle_ps(a1,a1,_MM_SHUFFLE( 1,0,3,2 ));
    b1 = _mm256_max_ps(a1,a2);
    a3 = _mm256_shuffle_ps(a1,a1,_MM_SHUFFLE( 0,3,2,1 ));
    a4 = _mm256_shuffle_ps(a1,a1,_MM_SHUFFLE( 2,1,0,3 ));
    b2 = _mm256_max_ps(a3,a4);
    c1 = _mm256_max_ps(b1,b2);
    c1= _mm256_max_ps(c1,_mm256_permute2f128_ps(c1,c1,_MM_SHUFFLE( 0,0,0,1 )));
    t = (float *) &c1;
	p.val = t[0];
	// determine max index
	unsigned int mask = _mm256_movemask_ps(_mm256_cmp_ps(ft,c1,_CMP_EQ_OQ));
	c1 = _mm256_cmp_ps(ft,c1,_CMP_EQ_OQ);
    p.idx = (int) (__builtin_ffs(mask));
    return p;
}

// Uses 64bit pdep / pext to save a step in unpacking.
// found on: https://stackoverflow.com/questions/36932240/avx2-what-is-the-most-efficient-way-to-pack-left-based-on-a-mask
#ifdef __AVX2__
__m256 compress256(__m256 src, unsigned int mask /* from movmskps */)
{ 
	uint64_t expanded_mask = _pdep_u64(mask, 0x0101010101010101);  // unpack each bit to a byte
	expanded_mask *= 0xFF;    // mask |= mask<<1 | mask<<2 | ... | mask<<7;
	// ABC... -> AAAAAAAABBBBBBBBCCCCCCCC...: replicate each bit to fill its byte

	const uint64_t identity_indices = 0x0706050403020100;    // the identity shuffle for vpermps, packed to one index per byte
	uint64_t wanted_indices = _pext_u64(identity_indices, expanded_mask);

	__m128i bytevec = _mm_cvtsi64_si128(wanted_indices);
	__m256i shufmask = _mm256_cvtepu8_epi32(bytevec);

	return _mm256_permutevar8x32_ps(src, shufmask);
}
#else
// slow and simple AVX implementation
__m256i compress256_AVX(__m256i vals, unsigned int mask) {
    int tmpIdx[8];
    memset(tmpIdx,0,8*sizeof(int));
    __m256i idxSelect;
    int *idx = (int*) &vals;
    int counter = 0;
    for (int k=0; k < 8; k++) {
        if ((mask & 1) != 0) {
            tmpIdx[counter++] = idx[k];
        }
        mask >>= 1;
    }
    idxSelect = _mm256_set_epi32(tmpIdx[7],tmpIdx[6],tmpIdx[5],tmpIdx[4],tmpIdx[3],tmpIdx[2],tmpIdx[1],tmpIdx[0]);
    return idxSelect;
}
#endif



// we want the exact max to limit the amount of added hyperplanes --> requires two runs for each grid point Y
void makeElementListExact(float* aGamma, float* bGamma, float* ftStore, float* X, int dim, int nH, int N, int* idxMax, int* elementListLocal, int* counterLocal) {
    int i, k, mask;
    __m256 ft, cmp, a, x, val1, max;
	__m256i idxSelect;	
	float ftMax;
	struct maxS p;
	__m256 idx = _mm256_set_ps(7,6,5,4,3,2,1,0);
	max = _mm256_set1_ps(-FLT_MAX);
    ftMax = -FLT_MAX;

	// find max
    for (i=0; i < nH-(nH%8); i+=8) {
        ft = _mm256_load_ps(bGamma + i);
        // ftTmp = bGamma[i] + aGamma[idxA]*X[idxB];
        for (k=0; k < dim; k++) {
   		    x = _mm256_set1_ps(*(X + N*k)); 
            a = _mm256_loadu_ps(aGamma + i + k*nH); 
#ifdef __AVX2__
			ft = _mm256_fmadd_ps(x,a,ft); // combined multiply+add
#else
            ft = _mm256_add_ps(ft,_mm256_mul_ps(x,a));
#endif			
        }
		// store ft
		_mm256_store_ps(ftStore+i,ft);
		max = _mm256_max_ps(ft,max);
	}
	// horizontal max
	p = avx_max_struct(max,max);
	// does not reflect the true Index --> dummy value --> change of index should be required 
	*idxMax = 0;
	// broadcast max back to avx register
	ftMax = p.val;
   	val1 = _mm256_set1_ps(ftMax-250);
	for (i=0; i < nH-(nH%8); i+=8) {
        // ftTmp > ftInnerMax
		ft = _mm256_load_ps(ftStore+i);
       	cmp =_mm256_cmp_ps(ft,val1,_CMP_GT_OQ);
		mask = _mm256_movemask_ps(cmp);
		
		// store in elementList
		if (mask != 0) {
#ifdef __AVX2__	
			idxSelect = _mm256_cvtps_epi32(compress256(_mm256_add_ps(idx,_mm256_set1_ps(i)), mask));
#else
			idxSelect = compress256_AVX(_mm256_cvtps_epi32(_mm256_add_ps(idx,_mm256_set1_ps(i))), mask);
#endif
			_mm256_storeu_si256((__m256i*) (elementListLocal+(*counterLocal)),idxSelect);
			int count = _mm_popcnt_u32(mask);
			*counterLocal += count;
		}
	}
	float ftInner;
	// for the remaining hyperplanes do it the scalar way
	for (i = i; i < nH; i++) {
    	ftInner = bGamma[i] + aGamma[i]**X;
		for (k=1; k < dim; k++) {
			ftInner += aGamma[i+k*nH]**(X + (k*N));
		}

		if (ftInner > ftMax-250) {
			/* find maximum element */
			if (ftInner > ftMax) {
				ftMax = ftInner;
				*idxMax = i;
			}

			elementListLocal[(*counterLocal)++] = i;
		}
	}
}

void makeElementListExactY(float* aGamma, float* bGamma, float* ftStore, float* X, int dim, int nH, int N, int* idxMax, int* elementListLocal, int* counterLocal, int* elementIds) {
    int i, k, mask;
    __m256 ft, cmp, a, x, val1, max;
    __m256i idxSelect;
    float ftMax;
    struct maxS p;
    max = _mm256_set1_ps(-FLT_MAX);
    ftMax = -FLT_MAX;

    // find max
    for (i=0; i < nH-(nH%8); i+=8) {
        ft = _mm256_load_ps(bGamma + i);
        // ftTmp = bGamma[i] + aGamma[idxA]*X[idxB];
        for (k=0; k < dim; k++) {
            x = _mm256_set1_ps(*(X + N*k));
            a = _mm256_loadu_ps(aGamma + i + k*nH);
#ifdef __AVX2__
            ft = _mm256_fmadd_ps(x,a,ft); // combined multiply+add
#else
            ft = _mm256_add_ps(ft,_mm256_mul_ps(x,a));
#endif
        }
        // store ft
        _mm256_store_ps(ftStore+i,ft);
        max = _mm256_max_ps(ft,max);
    }
    // horizontal max
    p = avx_max_struct(max,max);
    // does not reflect the true Index --> dummy value --> change of index should be required 
    *idxMax = 0;
    // broadcast max back to avx register
    ftMax = p.val;
    val1 = _mm256_set1_ps(ftMax-250);
    for (i=0; i < nH-(nH%8); i+=8) {
        // ftTmp > ftInnerMax
        ft = _mm256_load_ps(ftStore+i);
        cmp =_mm256_cmp_ps(ft,val1,_CMP_GT_OQ);
        mask = _mm256_movemask_ps(cmp);

        // store in elementList: elementListLocal[(*counterLocal)++] = elementIds[i]
        if (mask != 0) {
#ifdef __AVX2__ 
            idxSelect = _mm256_cvtps_epi32(compress256(_mm256_cvtepi32_ps(_mm256_loadu_si256((__m256i*) (elementIds+i))), mask));
#else
            idxSelect = compress256_AVX(_mm256_loadu_si256((__m256i*) (elementIds+i)), mask);
#endif
            _mm256_storeu_si256((__m256i*) (elementListLocal+(*counterLocal)),idxSelect);
            int count = _mm_popcnt_u32(mask);
            *counterLocal += count;
        }
    }
    float ftInner;
    // for the remaining hyperplanes do the scalar way
    for (i = i; i < nH; i++) {
        ftInner = bGamma[i] + aGamma[i]**X;
        for (k=1; k < dim; k++) {
            ftInner += aGamma[i+k*nH]**(X + (k*N));
        }

        if (ftInner > ftMax-250) {
            /* find maximum element */
            if (ftInner > ftMax) {
                ftMax = ftInner;
                *idxMax = i;
            }

            elementListLocal[(*counterLocal)++] = elementIds[i];
        }
    }
}



void inline findMaxVal(float* aGamma, float* bGamma, float* ftInner, float* X, int dim, int nH, int N, __m256* ftMax, float* boxEvalPoints, int* numPointsPerBox, int* idxElementsBox, int* numElementsBox) {
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

	// we include all hyperplanes in the list, that are less than 250 away from the max hyperplane (in log-scale)
    *ftMax = _mm256_add_ps(*ftMax,_mm256_set1_ps(-250));
    for (i=0; i < nH; i++) {
        ft = _mm256_load_ps(ftInner+8*i);
        for (k=0; k < dim; k++) {
            a = _mm256_set1_ps(*(aGamma + i*dim + k));
            evalTmp = _mm256_mul_ps(_mm256_sub_ps(a,aGammaMax[k]),sign[k]);
            cmp = _mm256_blendv_ps(zeros,ones,_mm256_cmp_ps(evalTmp,zeros,_CMP_GT_OQ));
            val1 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(evalTmp,Delta[k]),cmp),numPoints);
            ft = _mm256_add_ps(ft,val1);
        }
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


void preCondGradAVXC(int** elementList, int** elementListSize, int* numEntries, int* maxElement, int* idxEntries, float* X, float* grid, unsigned short int* YIdx, int *numPointsPerBox, float* boxEvalPoints, int numBoxes, double* a, double* aTrans, double* b, float gamma, float weight, float* delta, int N, int M, int dim, int nH) {
	/* initialize elementList */
    static const int elementListIncrement = 10000000;
	free(*elementListSize);
	free(*elementList);
    *elementListSize = malloc(sizeof(int)); **elementListSize = elementListIncrement;
	*elementList = malloc(**elementListSize*sizeof(int));

	float *aGamma, *aTransGamma, *bGamma;
	alloc_aligned_mem(dim*nH,ALIGN,&aGamma);
	alloc_aligned_mem(dim*nH,ALIGN,&aTransGamma);
	alloc_aligned_mem(nH,ALIGN,&bGamma);

	int i,j,k;
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
			aTransGamma[i+(k*nH)] = gamma*aTrans[i+(k*nH)];
		}
		bGamma[i] = gamma*b[i];
	}
    //double iStart = cpuSecond();

	counter = 0; savedValues = 0; savedValues2 = 0;
	#pragma omp parallel
	{   
		float* ftInner;
		alloc_aligned_mem(nH,ALIGN,&ftInner);
		int sizeElementList = elementListIncrement;
		int *elementListLocal = malloc(sizeElementList*sizeof(int));
		int *idxMax = malloc(sizeof(int));
		int *counterLocal = malloc(sizeof(int));
		int numElements, numElementsOld;
		/* calculate gradient for samples */
		*counterLocal = 0;
		#pragma omp for
		for (j=0; j < N; j++) {		
			if (sizeElementList < *counterLocal + nH) {
				sizeElementList *= 2;
				int *tmp = realloc(elementListLocal,sizeElementList*sizeof(int));
				if (tmp != NULL) {
				   elementListLocal = tmp; 
				} else {
					error("Realloc failed --> aborting execution");
				}
			}
			numElementsOld = *counterLocal;
			makeElementListExact(aGamma, bGamma, ftInner, X+j, dim, nH, N, idxMax, elementListLocal, counterLocal);
			numElements = *counterLocal-numElementsOld;
			maxElement[j] = *idxMax;
			numEntries[j] = numElements;
		}

		#pragma omp critical
		{   
			/* total number of elements that where added */
			counter += *counterLocal;
		}

		#pragma omp barrier

		/* reallocate idxList if neccessary */
		#pragma omp single
		{   
			if (**elementListSize < counter) {
				Rprintf("Reallocate elementList\n");
				
				int *tmp = realloc(*elementList,counter*sizeof(int));
				if (tmp != NULL) {
					*elementList = tmp;
				} else { 
					error("Realloc failed --> aborting execution");
				}
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
				memcpy((*elementList)+savedValues,elementListLocal,*counterLocal*sizeof(int));
				savedValues += *counterLocal;
			}
		}
		free(elementListLocal); free_aligned_mem(ftInner); free(idxMax); free(counterLocal);
	}
    //double timeTotal = cpuSecond()-iStart;

	double i1, part1, part2, part3, part4, part5;
	part1 = part2 = part3 = part4 = part5 = 0;

	//iStart = cpuSecond();
	#pragma omp parallel
	{
		float *stInner; 
		alloc_aligned_mem(8*nH,ALIGN,&stInner);

		int sizeElementList = elementListIncrement;
		int *elementListLocal = malloc(sizeElementList*sizeof(int));
        int *idxElements = malloc(8*nH*sizeof(int));
        int *idxElementsBox = malloc(8*nH*sizeof(int));
        int numElements, numElementsOld, numElementsBox[8], totalElements = 0;
		float *aLocal, *bLocal;
		alloc_aligned_mem(dim*nH,ALIGN,&aLocal);
		alloc_aligned_mem(nH,ALIGN,&bLocal);
		float *Ytmp = calloc(8*dim,sizeof(float));
		int *idxEntriesLocal = malloc(M*sizeof(int));
		int *maxElementLocal = malloc(M*sizeof(int));
		int *numEntriesLocal = malloc(M*sizeof(int));
	    int *idxMax = malloc(sizeof(int));
        int *counterLocal = malloc(sizeof(int));
		int l,n,m,p;
		__m256 stMax_;
		*counterLocal = 0;
		#pragma omp for schedule(dynamic,1) private(j,i,k,p)
        for (n = 0; n < numBoxes; n+=8) {
			/* check for active hyperplanes */
            /* eval all hyperplanes for some corner point of the box */
            i1 = cpuSecond();
            memset(Ytmp, 0, 8*dim*sizeof(float));
            for (m = 0; m < fmin(n+8,numBoxes)-n; m++) {
                for (k=0; k < dim; k++) {
                    Ytmp[m + k*8] = boxEvalPoints[(n+m)*3*dim + k];
                }
            }
            findMaxVal(aTransGamma, bGamma, stInner, Ytmp, dim, nH, 8, &stMax_,boxEvalPoints+n*3*dim,numPointsPerBox+n,idxElementsBox,numElementsBox);
            part1 += cpuSecond()-i1;
				
            for (m = 0; m < fmin(n+8,numBoxes)-n; m++) {
                j = n+m;
				//printf("%d,%d,%d (%d): %d nH, %d elements\n",n,m,j,omp_get_thread_num(),numElementsBox[m],numPointsPerBox[j+1] - numPointsPerBox[j]);
                i1 = cpuSecond();
                // save hyperplanes in one compact vector
                for (i=0; i < numElementsBox[m]; i++) {
                    bLocal[i] = bGamma[idxElementsBox[i+m*nH]];
                }
				for (k=0; k < dim; k++) {
                    for (i=0; i < numElementsBox[m]; i++) {
                        aLocal[i + k*numElementsBox[m]] = aGamma[idxElementsBox[i+m*nH]+k*nH];
                    }
                }
				part2 += cpuSecond() - i1;
	
				i1 = cpuSecond();
				int numPoints = numPointsPerBox[j+1]-numPointsPerBox[j];
				
				for (p=0; p < numPoints; p++) {
					l = p + numPointsPerBox[j];
					for (k=0; k < dim; k++) {
						Ytmp[k] = gridLocal[k]+delta[k]*YIdx[l*dim + k];
					}
					if (sizeElementList < *counterLocal + numElementsBox[m]) {
						sizeElementList *= 2;
						int *tmp = realloc(elementListLocal,sizeElementList*sizeof(int));
						if (tmp != NULL) {
						   elementListLocal = tmp; 
						} else {
							error("Realloc failed --> aborting execution");
						}
					}
					numElementsOld = *counterLocal;
					makeElementListExactY(aLocal, bLocal, stInner, Ytmp, dim, numElementsBox[m], 1, idxMax, elementListLocal, counterLocal, idxElementsBox+m*nH);
					numElements = *counterLocal-numElementsOld;
					
					maxElementLocal[totalElements] = idxElementsBox[*idxMax];
					numEntriesLocal[totalElements] = numElements;
					idxEntriesLocal[totalElements++] = l;
				}
				part3 += cpuSecond()-i1;
			}
		}
	    #pragma omp critical
        {   
            counter += *counterLocal; /* final number of elements in list */
        }
  		#pragma omp barrier
 
        // reallocate idxList if neccessary 
        #pragma omp single
        { 
         	if (**elementListSize < counter) {
   				int *tmp = realloc(*elementList,counter*sizeof(int));
				if (tmp != NULL) {
					*elementList = tmp;
				} else { 
					error("Realloc failed --> aborting execution");
				}
				**elementListSize = counter;
            }
        }

#ifdef _OPENMP
		int numThreads = omp_get_num_threads();
#else
		int numThreads = 1;
#endif

        // enforce ordered copying of memory
        #pragma omp for ordered schedule(static,1)
        for(j = 0; j < numThreads; j++)
        {
            #pragma omp ordered
            {
				memcpy((*elementList)+savedValues,elementListLocal,*counterLocal*sizeof(int));
				savedValues += *counterLocal;
				memcpy(idxEntries+savedValues2,idxEntriesLocal,totalElements*sizeof(int));
				memcpy(maxElement+savedValues2+N,maxElementLocal,totalElements*sizeof(int));
				memcpy(numEntries+savedValues2+N,numEntriesLocal,totalElements*sizeof(int));
				savedValues2 += totalElements;
			}
        }
  		free(Ytmp); free(stInner); free(idxElements); free(idxElementsBox); free(elementListLocal); free(idxEntriesLocal); free(maxElementLocal); free(numEntriesLocal); free(idxMax); free(counterLocal); free_aligned_mem(aLocal); free_aligned_mem(bLocal);

	} /* end of pragma parallel */
	**elementListSize = counter;
	free_aligned_mem(aGamma); free_aligned_mem(bGamma); free_aligned_mem(aTransGamma); 
}
#else
void preCondGradAVXC(int** elementList, int** elementListSize, int* numEntries, int* maxElement, int* idxEntries, float* X, float* grid, unsigned short int* YIdx, int *numPointsPerBox, float* boxEvalPoints, int numBoxes, double* a, double* aTrans, double* b, float gamma, float weight, float* delta, int N, int M, int dim, int nH) {
	// empty function if AVX is not used
}
#endif
