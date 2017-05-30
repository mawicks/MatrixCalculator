#include "constant.h"
#include <math.h>

void rotate ();
double twobytwo ();

#define sq(x) ((x)*(x))


void control1 (a, b, q, order, inputs, result)
double **a, **b, **q, *result;
int order, inputs;
{
  double v[2];
  double m11, m12, m22, norm;
  double m11max, m12max, m22max;
  double lastnorm = HUGE_VAL, tmp;
  int dropped = 1;
  int i, j, n;

/*
 * The distance to an uncontrollable pair
 * The real perturbation, rank 1 case.
 */

  while (dropped) {
    dropped = 0;
    for (i=0; i<order-1; i++) {
      m11 = m12 = m22 = 0.0;
      for (j=0; j<inputs; j++) {
	m11 += b[i][j] * b[i][j];
	m12 += b[i][j] * b[order-1][j];
	m22 += b[order-1][j] * b[order-1][j];
	addflops (6);
      }
      for (j=0; j<order-1; j++) {
	if (i!=j) {
	  m11 += a[i][j] * a[i][j];
	  m12 += a[i][j] * a[order-1][j];
	  m22 += a[order-1][j] * a[order-1][j];
	  addflops (6);
	}
      }
      m11 += sq (a[i][i] - a[order-1][order-1]);
      m12 += (a[i][i]-a[order-1][order-1])*a[order-1][i];
      m22 += sq (a[order-1][i]);
      m11 += sq (a[i][order-1]);
      addflops (10);
      if (m11 > m22) {
/*	v[0] = 0.5*(m11 - m22 + sqrt (sq(m11-m22)+4*sq(m12))); */
	v[0] = fabs (m11 - m22);
	v[1] = m12;
	tmp = sqrt (sq(v[0])+sq(v[1]));
	v[0] /= tmp; v[1] /= tmp;
/*	addflops (13); addsqrts (2); */
	addflops (6); addsqrts (2);
      }
      else {
	v[0] = 0.0; v[1] = 1.0;
      }
      rotate (a, b, q, v, order, inputs, i, order-1);
      if (m22 < lastnorm) {
	dropped = 1;
	lastnorm  = m22;
      }
    }
  }
  *result = sqrt (lastnorm);
  addsqrts (1);
}


void rotate (a, b, q, v, order, inputs, i, k)
double **a, **b, **q;
double v[2];
int order, inputs, i, k;
{
	int j;
	double tmp1, tmp2;

	/* Do reduction on b */
	for (j=0; j<inputs; j++) {
		tmp1 =  v[0]*b[i][j]+v[1]*b[k][j];
		tmp2 = -v[1]*b[i][j]+v[0]*b[k][j];
		b[i][j] = tmp1;
		b[k][j] = tmp2;
	}
	if (inputs > 0) addflops (inputs*6);
	/* Do reduction on a */
	for (j=0; j<order; j++) {
		tmp1 =  v[0]*a[i][j]+v[1]*a[k][j];
		tmp2 = -v[1]*a[i][j]+v[0]*a[k][j];
		a[i][j] = tmp1;
		a[k][j] = tmp2;
	}
	/* Multiply a again on right side */
	for (j=0; j<order; j++) {
		tmp1 =  v[0]*a[j][i]+v[1]*a[j][k];
		tmp2 = -v[1]*a[j][i]+v[0]*a[j][k];
		a[j][i] = tmp1;
		a[j][k] = tmp2;
	}
	if (order>0) addflops (order*6);
	/*
	 * Multiply q on right side to retain coordinate
	 * info
	 */
	for (j=0; j<order; j++) {
		tmp1 =  v[0]*q[j][i]+v[1]*q[j][k];
		tmp2 = -v[1]*q[j][i]+v[0]*q[j][k];
		q[j][i] = tmp1;
		q[j][k] = tmp2;
	}
	if (order>0) addflops (6*order);
}
