#include <R.h>

void printAVXInfo() {
#ifdef __AVX2__
	Rprintf("AVX2 and AXV vector extensions activated.\n");
#elif __AVX__
	Rprintf("AVX vector extensions activated.\n");
#else
	Rprintf("No vector extensions activated. (no AVX)\n");
#endif

#ifdef _OPENMP
	Rprintf("Parallel computing activated (OMP).\n");
#else
	Rprintf("No parallel computing activated (no OMP).\n");
#endif
}
