#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <sys/time.h>
#include <headers.h>
#ifdef _WIN32
#include <malloc.h>
#endif

/* Alligned memory allocation is non-trivial: We have seperate cases for Windows, Solaris and Linux/MacOS */
void alloc_aligned_mem(int numEntries, int align, float **ptr) {
#ifndef _WIN32
#ifdef __sun
	*ptr = memalign((size_t) align, (size_t) (numEntries*sizeof(float)));
#else
	int returnVal;
	returnVal = posix_memalign((void **) ptr, align, numEntries*sizeof(float));
	if (returnVal != 0) {
		error("Memory allocation failed\n");
	}
#endif
#else
	*ptr = _aligned_malloc(((size_t) numEntries)*sizeof(float),(size_t) align);
#endif
}

void free_aligned_mem(float *ptr) {
#ifndef _WIN32
	free(ptr);
#else
	_aligned_free(ptr);
#endif
}

void unzipParams(double *params, double *a, double *b, int dim, int nH, int transpose) {
    int i,j;
    if (transpose==1) {
        // transpose operation
        for (i=0; i < dim; i++) {
           for (j=0; j < nH; j++) {
               a[j*dim + i] = params[j + i*nH];
           }
        }
    } else {
        for (i=0; i < dim*nH; i++) {
            a[i] = params[i];
        }
    }
    for (i=0; i < nH; i++) {
        b[i] = params[dim*nH+i];
    }
}

void unzipParamsFloat(double *params, float *a, float *b, int dim, int nH, int transpose) {
    int i,j;
    if (transpose==1) {
        // transpose operation
        for (i=0; i < dim; i++) {
           for (j=0; j < nH; j++) {
               a[j*dim + i] = params[j + i*nH];
           }
        }
    } else {
        for (i=0; i < dim*nH; i++) {
            a[i] = params[i];
        }
    }
    for (i=0; i < nH; i++) {
        b[i] = params[dim*nH+i];
    }
}

double calcLambdaSq(double* grad, double* newtonStep, int dim, int nH) {
    double lambdaSq = 0;
    int i;
    for (i=0; i < nH*(dim+1); i++) {
        lambdaSq += grad[i]*-newtonStep[i];
    }
    return lambdaSq;
}

void copyVector(double* dest, double* source, int n, int switchSign) {
    int i;
    if (switchSign == 1) {
        for (i=0; i < n; i++) {
            dest[i] = -source[i];
        }
    } else {
        for (i=0; i < n; i++) {
            dest[i] = source[i];
        }
    }
}

double cpuSecond() {
    struct timeval tp;
    gettimeofday(&tp,NULL);
    return ((double)tp.tv_sec + (double)tp.tv_usec*1.e-6);
}


