#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <R_ext/Lapack.h>
#include <R.h>
#include <headers.h>

void recalcParamsC(double* X, double* y, int* T, int *lenT, int *dim_, double* aOptNew, double* bOptNew)
{
    int i;
	int dim = *dim_;
    #pragma omp parallel
    {   
        int j,k,info;
        double *A = malloc((dim+1)*(dim+1)*sizeof(double));
        double *b = malloc((dim+1)*sizeof(double));
        int dimT = dim+1;
        const int n = dim+1, nrhs=1, lda = dim+1;
        int *ipiv = malloc((dim+1)*sizeof(ptrdiff_t));
        #pragma omp for
        for (i=0; i < *lenT; i++) {
            for (j=0; j < dimT; j++) {
                b[j] = -y[T[i*dimT + j]];
            }

            for (j=0; j < dimT; j++) {
                for (k=0; k < dim; k++) {
                    A[j + dimT*k] = X[T[i*dimT + j]*dim + k];
                }
            }
            /* last row are ones */
            for (j=0; j < dimT; j++){
                A[dimT*dim+j] = 1;
            }

            F77_CALL(dgesv)(&n,&nrhs,A,&lda,ipiv,b,&n,&info);

            for (j=0; j < dim; j++) {

                aOptNew[i*dim + j] = b[j];
            }
            bOptNew[i] = b[dim];
        }
        free(A); free(b); free(ipiv);
    }
}

int comp(const void *a,const void *b)
{   
    double *x = (double *) a;
    double *y = (double *) b;

    if (*x < *y) return -1;
    else if (*x > *y) return 1; return 0;
}

double factorial(unsigned int f)
{   
    if ( f == 0 )
        return 1;
    return (double) (f * factorial(f - 1));
}


double calcGD(double* y, int lenY) {

    double eps;
    int i;

    eps = pow(10,-3);
    if (lenY == 1) {
        return exp(y[0]);
    } else if (y[lenY-1]-y[0] < eps) {
        /* make Taylor approximation */
        double yBar = 0;
        unsigned int d;
        double *v = malloc(lenY*sizeof(double));
/*      double v[lenY];*/
        double sumPow2 = 0; double sumPow3 = 0;

        for (i=0; i < lenY; i++) {
            yBar += y[i];
        }
        yBar = yBar/lenY;

        d = (unsigned int) lenY-1;
        for (i=0; i < lenY; i++) {
            v[i] = y[i]-yBar;
        }

        for (i=0; i < lenY; i++) {
            sumPow2 += v[i]*v[i];
            sumPow3 += v[i]*v[i]*v[i];
        }
        free(v);
        return exp(yBar)*(1/factorial(d) + 1/(2*factorial(d+2))*sumPow2 + 1/(3*factorial(d+3))*sumPow3);
    } else {
        double *y1 = malloc((lenY-1)*sizeof(double));
        double *y2 = malloc((lenY-1)*sizeof(double));
/*      double y1[lenY-1], y2[lenY-1];*/
        double returnVal;
        for (i=0; i < lenY-1; i++) {
           y1[i] = y[i];
           y2[i] = y[i+1];
        }
        qsort(y1,lenY-1,sizeof(double),comp);
        qsort(y2,lenY-1,sizeof(double),comp);

        returnVal = (calcGD(y2,lenY-1)-calcGD(y1,lenY-1))/(y[lenY-1]-y[0]);
        free(y1); free(y2);
        return returnVal;
    }
}


double calcIntegral(double* X, double* y, int lenY, int* T, int lenT, int dim, double* AdSave, double* GdSave) {
    double integral = 0;
    # pragma omp parallel 
    {   
        int i,j,k;
        double Gd, Ad;

        double integralLocal = 0;
        double *yTmp = malloc((dim+1)*sizeof(double));
        double *Xtmp = malloc(dim*dim*sizeof(double));
        int dimT = dim+1;
        double wkopt;
        double* work;
        int N = dim, n = dim, lda = dim, ldvl = dim, ldvr = dim, info, lwork;
        double *wr = malloc(N*sizeof(double));
        double *wi = malloc(N*sizeof(double));
        double *vl = malloc(N*dim*sizeof(double));
        double *vr = malloc(N*dim*sizeof(double));
        double wrCurr, wiCurr;

        #pragma omp for
        for (i=0; i < lenT; i++) {
            for (j=0; j < dim+1; j++) {
                yTmp[j] = y[T[i*dimT + j]];
            }
            qsort(yTmp,dim+1,sizeof(double),comp);
            Gd = calcGD(yTmp,dim+1);

            /* create temporary matrix (X(T(i,2:end)) - repmat(X(T(i,1),:),dim,1) */
            for (j=0; j < dim; j++) {
                for (k=0; k < dim; k++) {
                    Xtmp[j*dim + k] = X[T[i*dimT + j + 1]*dim + k] - X[T[i*dimT]*dim + k];
                }
            }
            /* Query and allocate the optimal workspace */
            lwork = -1;
            F77_CALL(dgeev)( "N", "N", &n, Xtmp, &lda, wr, wi, vl, &ldvl, vr, &ldvr, &wkopt, &lwork, &info );
            lwork = (int)wkopt;
            work = (double*)malloc( lwork*sizeof(double) );
            /* Solve eigenproblem */
            F77_CALL(dgeev)( "N", "N", &n, Xtmp, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, &info );

            /* calculate product of (potentially) complex eigenvalues */
            for (j=1; j < dim; j++) {
                wrCurr = wr[j]; wiCurr = wi[j];
                wr[j] = wr[j-1]*wrCurr - wi[j-1]*wiCurr;
                wi[j] = wr[j-1]*wiCurr + wi[j-1]*wrCurr;
            }
            Ad = fabs(wr[dim-1]);
            integralLocal += Ad*Gd;
            free(work);
            /* the volume can be calculated from Ad by multiplying the factor 1/dim! */
            AdSave[i] = Ad; GdSave[i] = Gd;
        }
        #pragma omp critical
        {   
            integral += integralLocal;
        }

        free(yTmp); free(Xtmp); free(wr); free(wi); free(vl); free(vr);
    }
    return integral;
}


void calcExactIntegralC(double *X, double *y, int *T, int *lenT, int *lenY, int *dim, double *targetIntegral, double *intEps, double *a, double *b) {

    int i;
    int maxIter = 50;
    int iter = 0;
    double factor;
    double *yTmp = malloc((*lenY)*sizeof(double));
    double oldDistAbs, oldDist;
    int cont;
	double *Ad = malloc(*lenT*sizeof(double));
	double *Gd = malloc(*lenT*sizeof(double));
	double integral;

    integral = calcIntegral(X,y,*lenY,T,*lenT,*dim,Ad,Gd);
    while (fabs(*targetIntegral - integral) > *intEps && iter < maxIter) {
        factor = 1; cont = 1;
        oldDistAbs = fabs(*targetIntegral - integral);
        oldDist = *targetIntegral - integral;

        while (cont==1) {
            for (i=0; i < *lenY; i++) {
                yTmp[i] = y[i] + oldDist*factor;
            }
            //printf("y: %.4f\n",yTmp[0]);

            integral = calcIntegral(X,yTmp,*lenY,T,*lenT,*dim,Ad,Gd);
            if (fabs(*targetIntegral - integral) < oldDistAbs)
            {
                for (i=0; i < *lenY; i++) {
                    y[i] += oldDist*factor;
                }
                cont = 0;
            } else {
                //printf("\t NewInt: %.6f,%.3f\n",integral,factor);
                factor = factor/2;
                iter++;
                if (iter > 10) {
                    break;
                }
            }
        }
    }
	recalcParamsC(X,y,T,lenT,dim,a,b);

    free(yTmp); free(Gd); free(Ad);
}

