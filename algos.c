#include "constant.h"
#include "functions.h"
#include <math.h>
#include <stdio.h>

double algo1 (q, b, a, order, inputs, norm)
double **q, **b, **a, *norm;
{
	double **ab, **u, **uta, **utb;
	int i,j, done;
	ab = array (order, order+inputs);
	u = array (order, order);
	uta = array (order, order);
	utb = array (order, inputs);
	*norm = HUGE_VAL;
	do {
	  done = 1;
  	  for (i=0; i<order; i++) {
	    for (j=0; j<order; j++) {
	      ab[i][j] = a[i][j];
	      if (i == j) ab[i][j] -= a[order-1][order-1];
	    }
	    for (j=0; j<inputs; j++) 
	      ab[i][j+order] = b[i][j];
	  }
	  if (order > 0) addflops (order);
	  m_svd (u, ab, (double *) 0, order, order+inputs);
	  tm_mul (u, a, uta, order, order, order);
	  mm_mul (uta, u, a, order, order, order);

	  tm_mul (u, b, utb, order, order, inputs);
	  m_copy (utb, b, order, inputs);

	  mm_mul (q, u, uta, order, order, order);
	  m_copy (uta, q, order, order, order);
	  if (ab[order-1][order-1] < *norm) {
		done = 0;
		*norm = ab[order-1][order-1];
	  }
	  printf ("Norm = %g  lambda = %g\n", *norm, a[order-1][order-1]);
	}
	while (!done);
	free_array (ab);
	free_array (uta);
	free_array (utb);
	free_array (u);
}

double algo2 (q, b, a, order, inputs, norm)
double **q, **b, **a, *norm;
{
	double **ab, **u, **uta, **utb;
	int i,j, done;
	ab = array (order, order+inputs);
	u = array (order, order);
	uta = array (order, order);
	utb = array (order, inputs);
	*norm = HUGE_VAL;
	do {
	  done = 1;
  	  for (i=0; i<order; i++) {
	    for (j=0; j<order; j++) {
	      ab[i][j] = a[i][j];
	      if (i == j) ab[i][j] -= a[order-1][order-1];
	    }
	    for (j=0; j<inputs; j++) 
	      ab[i][j+order] = b[i][j];
	  }
	  if (order > 0) addflops (order);
	  m_lq (ab, (double *) 0, order, order+inputs);
	  m_qr (ab, u, order, order+inputs);
	  tm_mul (u, a, uta, order, order, order);
	  mm_mul (uta, u, a, order, order, order);

	  tm_mul (u, b, utb, order, order, inputs);
	  m_copy (utb, b, order, inputs);

	  mm_mul (q, u, uta, order, order, order);
	  m_copy (uta, q, order, order, order);
	  if (ab[order-1][order-1] < *norm) {
		done = 0;
		*norm = ab[order-1][order-1];
	  }
	  printf ("Norm = %g  lambda = %g\n", *norm, a[order-1][order-1]);
	}
	while (!done);
	free_array (ab);
	free_array (uta);
	free_array (utb);
	free_array (u);
}
			
