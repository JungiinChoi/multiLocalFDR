#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#if defined(__AVX__)  // x86 AVX
#include <immintrin.h>
#elif defined(__ARM_NEON__)  // ARM NEON
#include <arm_neon.h>
#endif

void CNS(double* s_k, double* y_k, double* sy, double* syInv, double step, double* grad, double* gradOld, 
         double* newtonStep, int numIter, int activeCol, int nH, int m) {
  
  double normTmp = 0, dotProd = 0, dotProd2 = 0;
  double s_k_tmp, y_k_tmp, t = 0, H0;
  double* gammaBFGS = (double*) malloc(nH * sizeof(double));
  double* q = (double*) malloc(nH * sizeof(double));
  double alphaBFGS[numIter];
  double betaBFGS, tmp;
  int iterVec[numIter];
  int activeIdx = activeCol * nH;
  
  int j, l, curCol;
  
  // Compute gradients
#pragma omp parallel for reduction(+:normTmp, dotProd, t) private(s_k_tmp)
  for (j = 0; j < nH; j++) {
    gammaBFGS[j] = grad[j] - gradOld[j];
    s_k_tmp = step * newtonStep[j];
    dotProd += gammaBFGS[j] * s_k_tmp;
    normTmp += s_k_tmp * s_k_tmp;
    t += gradOld[j] * gradOld[j];
    s_k[activeIdx + j] = s_k_tmp;
  }
  
  t = sqrtf(t);
  if (-dotProd / normTmp > 0) {
    t += -dotProd / normTmp;
  }
  
  dotProd = dotProd2 = 0;
  for (j = 0; j < nH; j++) {
    y_k_tmp = gammaBFGS[j] + t * s_k[activeIdx + j];
    y_k[activeIdx + j] = y_k_tmp;
    dotProd += y_k_tmp * s_k[activeIdx + j];
    dotProd2 += y_k_tmp * y_k_tmp;
  }
  
  sy[activeCol] = dotProd;
  syInv[activeCol] = 1 / sy[activeCol];
  H0 = sy[activeCol] / dotProd2;
  
  for (j = 0; j < numIter; j++) {
    iterVec[j] = activeCol - j;
    if (iterVec[j] < 0) {
      iterVec[j] = m + iterVec[j];
    }
  }
  
  // Initialize q with grad
  memcpy(q, grad, nH * sizeof(double));
  
  // First loop: Apply updates
  for (l = 0; l < numIter; l++) {
    curCol = iterVec[l];
    tmp = 0;
    
#pragma omp parallel for reduction(+:tmp)
    for (j = 0; j < nH; j++) {
      tmp += s_k[curCol * nH + j] * q[j];
    }
    
    alphaBFGS[curCol] = tmp * syInv[curCol];
    
#if defined(__AVX__)
    __m256d alpha = _mm256_set1_pd(alphaBFGS[curCol]);
#pragma omp parallel for
    for (j = 0; j < nH - nH % 4; j += 4) {
      __m256d q_ = _mm256_loadu_pd(q + j);
      q_ = _mm256_sub_pd(q_, _mm256_mul_pd(alpha, _mm256_loadu_pd(y_k + curCol * nH + j)));
      _mm256_storeu_pd(q + j, q_);
    }
    for (j = nH - nH % 4; j < nH; j++) {
      q[j] -= alphaBFGS[curCol] * y_k[curCol * nH + j];
    }
    
#elif defined(__ARM_NEON__)
    float64x2_t alpha_neon = vdupq_n_f64(alphaBFGS[curCol]);
#pragma omp parallel for
    for (j = 0; j < nH - nH % 2; j += 2) {
      float64x2_t q_ = vld1q_f64(q + j);
      q_ = vsubq_f64(q_, vmulq_f64(alpha_neon, vld1q_f64(y_k + curCol * nH + j)));
      vst1q_f64(q + j, q_);
    }
    for (j = nH - nH % 2; j < nH; j++) {
      q[j] -= alphaBFGS[curCol] * y_k[curCol * nH + j];
    }
    
#else
#pragma omp parallel for
    for (j = 0; j < nH; j++) {
      q[j] -= alphaBFGS[curCol] * y_k[curCol * nH + j];
    }
#endif
  }
  
  // Apply H0 scaling
  for (j = 0; j < nH; j++) {
    q[j] = H0 * q[j];
  }
  
  // Second loop: Final adjustment
  for (l = 0; l < numIter; l++) {
    curCol = iterVec[numIter - 1 - l];
    betaBFGS = 0;
    
#pragma omp parallel for reduction(+:betaBFGS)
    for (j = 0; j < nH; j++) {
      betaBFGS += y_k[curCol * nH + j] * q[j];
    }
    betaBFGS *= syInv[curCol];
    tmp = (alphaBFGS[curCol] - betaBFGS);
    
#if defined(__AVX__)
    __m256d tmp_ = _mm256_set1_pd(tmp);
#pragma omp parallel for
    for (j = 0; j < nH - nH % 4; j += 4) {
      __m256d q_ = _mm256_loadu_pd(q + j);
      q_ = _mm256_add_pd(q_, _mm256_mul_pd(tmp_, _mm256_loadu_pd(s_k + curCol * nH + j)));
      _mm256_storeu_pd(q + j, q_);
    }
    for (j = nH - nH % 4; j < nH; j++) {
      q[j] += s_k[curCol * nH + j] * tmp;
    }
    
#elif defined(__ARM_NEON__)
    float64x2_t tmp_neon = vdupq_n_f64(tmp);
#pragma omp parallel for
    for (j = 0; j < nH - nH % 2; j += 2) {
      float64x2_t q_ = vld1q_f64(q + j);
      q_ = vaddq_f64(q_, vmulq_f64(tmp_neon, vld1q_f64(s_k + curCol * nH + j)));
      vst1q_f64(q + j, q_);
    }
    for (j = nH - nH % 2; j < nH; j++) {
      q[j] += s_k[curCol * nH + j] * tmp;
    }
    
#else
#pragma omp parallel for
    for (j = 0; j < nH; j++) {
      q[j] += s_k[curCol * nH + j] * tmp;
    }
#endif
  }
  
  // Update Newton step
  for (j = 0; j < nH; j++) {
    newtonStep[j] = -q[j];
  }
  
  free(gammaBFGS);
  free(q);
}
