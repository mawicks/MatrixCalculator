#include "constant.h"
#include <math.h>

void rotate ();
double twobytwo ();

#define sq(x) ((x)*(x))


void control2 (a, b, q, order, inputs, result)
double **a, **b, **q, *result;
int order, inputs;
{
  double v[2];
  double m11, m12, m22, norm;
  double m11max, m12max, m22max;
  double lastnorm = HUGE_VAL, tmp;
  int dropped = 1;
  int i, j, imax, nmax, n, n1, n2;

/*
 * The distance to an uncontrollable pair
 * The real perturbation, rank 2 case.
 */

  while (dropped) {
    dropped = 0;
    for (i=0; i<order-2; i++) {
      for (n=1; n<=2; n++) {
	if ( n == 1) {
	  n1 = order - 1;
	  n2 = order - 2;
	}
	else {
	  n1 = order - 2;
	  n2 = order - 1;
	}
	m11 = m12 = m22 = 0.0;
	for (j=0; j<inputs; j++) {
	  m11 += b[i][j] * b[i][j];

	  m12 += b[i][j] * b[n1][j];
	  m22 += b[n1][j] * b[n1][j];
	  addflops (6);
	}
	for (j=0; j<order-2; j++) {
	  if (i!=j) {
	    m11 += a[i][j] * a[i][j];
	    m12 += a[i][j] * a[n1][j];
	    m22 += a[n1][j] * a[n1][j];
	    addflops (6);
	  }
	}
	m11 += sq(a[n2][n1]);
	m12 += - a[n2][i] * a[n2][n1];
	m22 += sq(a[n2][i]);
	m11 += sq (a[i][i] - a[n1][n1]);
	m12 += (a[i][i]-a[n1][n1])*a[n1][i];
	m22 += sq (a[n1][i]);
	m11 += sq (a[i][n1]);
	addflops (16);
	v[0] = 0.5*(m11 - m22 + sqrt (sq(m11-m22)+4*sq(m12)));
	v[1] = m12;
	tmp = sqrt (sq(v[0])+sq(v[1]));
	v[0] /= tmp; v[1] /= tmp;
	addflops (13); addsqrts (2);
	rotate (a, b, q, v, order, inputs, i, n1);
	/* Compute norm */
	norm = 0.0;
	for (j=0; j<inputs; j++) {
	  norm += sq(b[order-1][j]) + sq(b[order-2][j]);
	  addflops (4);
	}
	for (j=0; j<order-2; j++) {
	  norm += sq(a[order-1][j]) + sq(a[order-2][j]);
	  addflops (4);
	}
	if (norm < lastnorm) {
	  dropped = 1;
	  lastnorm  = norm;
	}
      }
    }
  }
  *result = sqrt (lastnorm);
  addsqrts (1);
}
