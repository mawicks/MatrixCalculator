#include <math.h>
#include "constant.h"
#include "functions.h"

void eigv (a, q, size)
double **a, **q;
int size;
{
	int is, il, first=0, last, shortest, longest = size;
	ident (q, size, size);
	m_hes (a, q, size);
	while (longest > 1) {
	last=size-1, shortest=size-1, longest = 0;
	shortest = size-1;
	il = 0; is = 0;
	while (is < size-1) {
		while (is<size-1 && 
		       (fabs(a[is+1][is])>EPS*fabs(a[is][is]) ||
		       fabs(a[is+1][is])>EPS*fabs(a[is+1][is+1])))
			is++;  /* "is" is next zero sub-diag element */
			addflops (1);
		if (is<size-1) 
			a[is+1][is] = 0.0; /* Prevent further pivots */
		if ((is-il > 1) && (is-il) < shortest) {
			first = il;
			last = is;
			shortest = is - il;
		}
		if (is-il > longest) longest = is - il;
		il = ++is;
	}
	if (longest > 1) {
		m_shft (a, q, size, first, last);
		m_hes (a, q, size);
	}
	}
	return;
}


		
			
			
