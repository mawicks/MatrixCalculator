#include "mem.h"
#include <stdlib.h>

double **array (m, n)
int m,n;
{
	double **vecs, **index;
	index = vecs = MALLOC (m, double *);
	*index = MALLOC (m*n, double);
	for (m--; m>0; m--, index++)
		*(index+1) = (*index) + n;
	return (vecs);
}

void free_array (a)
double **a;
{
	free (*a);
	free(a);
}
