#include <math.h>
#include <stdlib.h>
#include <string.h>

void calcKernelDens(double *X, double *sampleWeights, double *kernelDens, double *h, int *N_, int *d_) {

	int i,N,d;
	N = *N_;
	d = *d_;
	const double pi = 3.141592653589793115997963468544185161590576171875;
	double *weights = malloc(d*sizeof(double));
	double *hInv = malloc(d*sizeof(double));
	double normalization = 1;
	/* calculate normalization constants */
	for (i=0; i < d; i++) {
		weights[i] = 1/(h[i]*sqrt(2*pi));
		normalization *= weights[i]; 
		hInv[i] = -0.5/(h[i]*h[i]);
	}

	#pragma omp parallel
	{   
		int idxX,idxY,j,k;
		double innerSum,yTmp;
		/* for each point in Y calculate kernel density */
		#pragma omp for
		for (i=0; i < N; i++) {
			idxX = i*d;
			/* evalute the Gauss kernel for each point in X and multiply with sampleWeight */
			yTmp = 0;
			for (j = 0; j < N; j++) {
				idxY = j*d;
				/* for each point we need to evaluate 1-D Gauss Kernels and multiply them */
				innerSum = 0;
				for (k = 0; k < d; k++) {
					innerSum += (X[idxX+k]-X[idxY+k])*(X[idxX+k]-X[idxY+k])*hInv[k];
				}
				yTmp += sampleWeights[j]*exp(innerSum)*normalization;
			}
			kernelDens[i] = yTmp;
		}
	}

	free(weights); free(hInv);
}
