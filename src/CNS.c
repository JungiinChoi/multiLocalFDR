#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <immintrin.h>
#include <stdio.h>

void CNS(double* s_k, double *y_k, double *sy, double *syInv, double step, double *grad, double *gradOld, double *newtonStep, int numIter, int activeCol, int nH, int m) {

	double normTmp, dotProd, dotProd2;
	double s_k_tmp, y_k_tmp, t, H0;
	double* gammaBFGS = (double *) malloc(nH*sizeof(double));
	double* q = (double *) malloc(nH*sizeof(double));
	double alphaBFGS[numIter];
	double betaBFGS,tmp;
	int iterVec[numIter];
	int activeIdx = activeCol*nH; 	
	int j,l,curCol;

	normTmp = dotProd = t = 0;
	
	#pragma omp parallel for reduction(+:normTmp,dotProd,t) private(s_k_tmp)
	for (j=0; j < nH; j++) {
		gammaBFGS[j] = grad[j] - gradOld[j];
		s_k_tmp = step*newtonStep[j];
		dotProd += gammaBFGS[j]*s_k_tmp;
		normTmp += s_k_tmp*s_k_tmp;
		t += gradOld[j]*gradOld[j];
		s_k[activeIdx + j] = s_k_tmp;
	}

	t = sqrtf(t); // finish calculation of norm
	if (-dotProd/normTmp > 0) {
		t += -dotProd/normTmp;
	}
	//printf("gammaBFGS: %.4e, t: %.7f, gammaBFGS'*s_k: %.4e, s_k'*s_k: %.4e\n",gammaBFGS[0],t,dotProd,normTmp);

	dotProd = dotProd2 = 0;
	#pragma omp parallel for reduction(+:dotProd,dotProd2) private(y_k_tmp)
	for (j=0; j < nH; j++) {
		y_k_tmp = gammaBFGS[j] + t*s_k[activeIdx+j];
		y_k[activeIdx + j] = y_k_tmp;
		dotProd += y_k_tmp*s_k[activeIdx+j];
		dotProd2 += y_k_tmp*y_k_tmp;
	}
	sy[activeCol] = dotProd; 
	syInv[activeCol] = 1/sy[activeCol];

    H0 = sy[activeCol]/dotProd2;

	for (j=0; j < numIter; j++) {
		iterVec[j] = activeCol-j;
		if (iterVec[j] < 0) {
			iterVec[j] = m + iterVec[j];
		}
	}
	// q = grad
	memcpy(q,grad,nH*sizeof(double));
#ifdef __AVX__
	__m256d alpha,q_;
#endif
    // first for-loop
	for (l=0; l < numIter; l++) {
		curCol = iterVec[l];
		tmp = 0;
		#pragma omp parallel for reduction(+:tmp)
		//for (j=nH-nH%4; j < nH; j++) {
		for (j=0; j < nH; j++) {
			tmp += s_k[curCol*nH + j]*q[j];
		}
		alphaBFGS[curCol] = tmp*syInv[curCol];
#ifdef __AVX__
		alpha = _mm256_set1_pd(alphaBFGS[curCol]);
		#pragma omp parallel for private(q_)
		for (j=0; j < nH-nH%4; j+=4) {
			q_ = _mm256_loadu_pd(q+j);
			q_ = _mm256_sub_pd(q_,_mm256_mul_pd(alpha,_mm256_loadu_pd(y_k+curCol*nH+j)));
			_mm256_storeu_pd(q+j,q_);
		}
		for (j=nH-nH%4; j < nH; j+=4) {
			q[j] -= alphaBFGS[curCol]*y_k[curCol*nH+j];
		}
#else
		#pragma omp parallel for
		for (j=0; j < nH; j++) {
			q[j] -= alphaBFGS[curCol]*y_k[curCol*nH+j];
		}
#endif
	}

	for (j=0; j < nH; j++) {
		q[j] = H0*q[j]; // is "r" in the matlab code and the book
	}
	// second for-loop
#ifdef __AVX__
	__m256d tmp_;
#endif
	for (l=0; l < numIter; l++) {
		curCol = iterVec[numIter-1-l];
		betaBFGS = 0;
		#pragma omp parallel for reduction(+:betaBFGS)
		for (j=0; j < nH; j++) {
			betaBFGS += y_k[curCol*nH+j]*q[j];
		}
		betaBFGS *= syInv[curCol];
		tmp = (alphaBFGS[curCol]-betaBFGS);
#ifdef __AVX__
		tmp_ = _mm256_set1_pd(tmp);
		#pragma omp parallel for private(q_)
		for (j=0; j < nH-nH%4; j+=4) {
			q_ = _mm256_loadu_pd(q+j);
			q_ = _mm256_add_pd(q_,_mm256_mul_pd(tmp_,_mm256_loadu_pd(s_k+curCol*nH+j)));
			_mm256_storeu_pd(q+j,q_);
		}
		for (j=nH-nH%4; j < nH; j++) {
			q[j] += s_k[curCol*nH + j]*tmp;
		}
#else
		#pragma omp parallel for
		for (j=0; j < nH; j++) {
			q[j] += s_k[curCol*nH + j]*tmp;
		}
#endif
	}

	for (j=0; j < nH; j++) {
		newtonStep[j] = -q[j];
	}
	free(gammaBFGS); free(q);
}
