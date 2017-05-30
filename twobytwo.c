#include <math.h>
#include "constant.h"
#define sq(x) ((x)*(x))

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
   
#define fixup(r,s,mx)					\
	     ( (s != 0.0) ?				\
	       ( mx = max(fabs(r),fabs(s)),		\
		 r /= mx,				\
                 s /= mx,				\
	         mx = sqrt (r*r + s*s),			\
 	         r /= mx,		 	 	\
		 s /= mx,                               \
		 addflops(7),                           \
		 addsqrts(1),                           \
		 1):			        	\
               0                                        \
 	      )

double twobytwo (a, b, v)
double a[2][2], b[2][2], v[2];
{
  /* This subroutine solves the distance to an uncontrollable pair
     problem for the two by two matrices, A and B.  The general 
     problem, for both real and complex perturbations, can be
     reduced to applications of this subproblem to pairs of rows
     of the more general matrices. */

  double last_norm = HUGE_VAL, r, s, tmp1, tmp2;
  double m[3];
  int i, count=0;
  v[0] = 1.0;
  v[1] = 0.0;

  while ((count < 3) &&
         (m[2] = sq(b[1][0]) + sq(b[1][1]) + sq(a[1][0])) < last_norm) {
    addflops (5);
    count += 1;
    last_norm = m[2];
    m[0] = sq(b[0][0]) + sq(b[0][1]) + sq(a[0][1]);
    addflops (5);
    if (m[0] >= m[2]) {
	m[0] += sq(a[0][0]-a[1][1]);
    	m[1] = b[0][0]*b[1][0] + b[0][1]*b[1][1] + (a[0][0]-a[1][1])*a[1][0];
	r = fabs (m[0] - m[2]);
   	s = m[1];
	addflops (10);
    }
    else {
        r = 0.0;
	s = 1.0;
    }
    if (!fixup (r, s, tmp1)) break;
    for (i=0; i<2; i++) {	/* Multiply B on the left  */
      tmp1 =  r*b[0][i] + s*b[1][i];
      tmp2 = -s*b[0][i] + r*b[1][i];
	addflops (6);
      b[0][i] = tmp1;
      b[1][i] = tmp2;
    }
    for (i=0; i<2; i++) {	/* Multiply A on the left  */
      tmp1 =  r*a[0][i] + s*a[1][i];
      tmp2 = -s*a[0][i] + r*a[1][i];
       addflops (6);
      a[0][i] = tmp1;
      a[1][i] = tmp2;
    }
    for (i=0; i<2; i++) {	/* Multiply A on the right */
      tmp1 =  r*a[i][0] + s*a[i][1];
      tmp2 = -s*a[i][0] + r*a[i][1];
	addflops (6);
      a[i][0] = tmp1;
      a[i][1] = tmp2;
    }
				/* Multply V on the left */
      tmp1 =  r*v[0] - s*v[1];
      tmp2 = -s*v[0] - r*v[1];
	addflops (6);
      v[0] =   tmp1;
      v[1] = - tmp2;
  }
  return (last_norm);
}



