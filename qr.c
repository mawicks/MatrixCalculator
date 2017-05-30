#include <math.h>
#include <stdio.h>

#include "constant.h"
#include "functions.h"
#include "mem.h"
#define min(a,b) ((a)<(b) ? (a) : (b))
#define max(a,b) ((a)>(b) ? (a) : (b))
#define swap(a,b,tmp) (tmp=a,a=b,b=tmp)

/* The following are two ugly macros which simplify
   coding while maintaining execution speed          */
   
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

#define rotate(a1,a2,r,s,tmp)			\
            ( tmp =  r * a1 + s * a2,		\
	      a2  = -s * a1 + r * a2,		\
              a1  = tmp,                        \
	      addflops(6)                       \
	     )

void m_qr (a, q, rows, cols)
/************************************************************************
 *  Purpose:  This subroutine computes a QR factorization of a 	 	*
 *            matrix using orthogonal Givens rotations			*
 * 									*
 *  Inputs:   a     - (double) array containing matrix to be factored	*
 *            rows  - (integer) number of rows in a		     	*
 *            cols  - (integer) number of columns in a		     	*
 *  Outputs:  a     - (double) upper triangular matrix replaces a	*
 *            q     - (double) orthogonal matrix overwrites q		*
 ************************************************************************
 */									
double **a, **q;
int rows, cols;
{
   int i, j, k, clim;
   double r, s, tmp;
   
   if (q != 0) ident (q, rows, rows);   
   clim = min (rows, cols);
   for (i=0; i<clim; i++) {		/*  For each column */
      for (j=i+1; j<rows; j++) {	/*  Go down the column */
         r = a[i][i];
	 s = a[j][i];
	 if (fixup (r, s, tmp)) {
	    for (k=i; k<cols; k++) { 	/*  Do elimination */
	       rotate (a[i][k],a[j][k], r, s, tmp);
	    }                           /*  Closes k loop */ 
	    if (q != 0) for (k=0; k<rows; k++) { /*  Do multiplication on Q */
	       rotate (q[k][i],q[k][j], r, s, tmp);
	    }
	    a[j][i] = 0.0;  		/*  Kill it */
         }				/*  Closes if */
      }         			/*  Closes j loop */
   }					/*  Closes i loop */
}

void m_lq (a, q, rows, cols)
/************************************************************************
 *  Purpose:  This subroutine computes an LQ factorization of a	 	*
 *            matrix using orthogonal Givens rotations			*
 * 									*
 *  Inputs:   a     - (double) array containing matrix to be factored	*
 *            rows  - (integer) number of rows in a		     	*
 *            cols  - (integer) number of columns in a		     	*
 *  Outputs:  a     - (double) upper triangular matrix replaces a	*
 *            q     - (double) orthogonal matrix overwrites q		*
 ************************************************************************
 */									
double **a, **q;
int rows, cols;
{
   int i, j, k, clim;
   double r, s, tmp;
   
   if (q != 0) ident (q, rows, rows);   
   clim = min (rows, cols);
   for (i=0; i<clim; i++) {		/*  For each row */
      for (j=i+1; j<cols; j++) {	/*  Go across the column */
         r = a[i][i];
	 s = a[i][j];
         if (fixup (r, s, tmp)) {
	    for (k=i; k<rows; k++) { 	/*  Do elimination */
	       rotate (a[k][i], a[k][j], r, s, tmp);
	    }                           /*  Closes k loop */ 
	    if (q != 0) for (k=0; k<cols; k++) { /*  Do multiplication on Q */
	       rotate (q[i][k], q[j][k], r, s, tmp);
	    }
	    a[i][j] = 0.0;  		/*  Kill it */
         }				/*  Closes if */
      }         			/*  Closes j loop */
   }					/*  Closes i loop */
}


void m_bidiag (q, a, v, rows, cols)

/************************************************************************
 *  Purpose:  This subroutine bidiagonalizes a matrix   	 	*
 *            using orthogonal Givens rotations		     		*
 *            The matrix is factored as q * bidiag * v                  *
 * 									*
 *  Inputs:   a     - (double) array containing matrix to be factored	*
 *            rows  - (integer) number of rows in a		     	*
 *            cols  - (integer) number of columns in a		     	*
 *  Outputs:  a     - (double) upper triangular matrix replaces a	*
 *            q     - (double) orthogonal matrix overwrites q		*
 *            v     - (double) orthogonal matrix overwrites v           *
 ************************************************************************
 */									
double **a, **q, **v;
int rows, cols;
{
   int i, ipiv, j, k, clim;
   double r, s, tmp;
   
   if (q != 0) ident (q, rows, rows);   
   if (v != 0) ident (v, cols, cols);   
   if (cols > rows) 
      m_lq (a, v, rows, cols); 		/*  Necessary to reduce size */
   clim = min (rows, cols);
   for (i=0; i<clim; i++) {		/*  For each diagonal element */
      for (j=i+1; j<rows; j++) {	/*  Go down the column */
         ipiv = i;
         r = a[ipiv][i];
	 s = a[j][i];
	 if (fixup (r, s, tmp)) {
	    for (k=i; k<clim; k++) { 	/*  Do elimination */
	       rotate (a[ipiv][k], a[j][k], r, s, tmp);
	    }                           /*  Closes k loop */ 
	    if (q != 0) for (k=0; k<rows; k++) { /*  Do multiplication on Q */
	       rotate (q[k][ipiv], q[k][j], r, s, tmp);
	    }
	    a[j][i] = 0.0;  		/*  Kill it */
         }				/*  Closes if */
      }         			/*  Closes j loop */
      for (j=i+2; j<clim; j++) {	/*  Go across the row */
         ipiv = i+1;
         r = a[i][ipiv];
	 s = a[i][j];
	 if (fixup (r, s, tmp)) {
	    for (k=i; k<rows; k++) {    /*  Multiply columns */
	       rotate (a[k][ipiv], a[k][j], r, s, tmp);
	    }      
	    if (v != 0) for (k=0; k<cols; k++) { /*  Multiply rows of V */
	       rotate (v[ipiv][k], v[j][k], r, s, tmp);
	    }
	    a[i][j] = 0.0;		/*  Kill it */	    
	 }				/*  Close if  */
      }					/*  Closes j loop */
   }					/*  Closes i loop */      
}

void m_svd (u, a, v, rows, cols)
double **a, **v, **u;
int rows, cols;
{
   double tmp, r, s;
   double u1, u2, u3, l2;
   int i, k, lim, start = 0, end;
   lim = min (rows, cols);
   end = lim - 2;
   m_bidiag (u, a, v, rows, cols);
   while (start < lim-1) {
/* Find beginning of first block */
     while (start < lim-1 &&
            fabs(a[start][start+1]) <=
	    EPS*(fabs(a[start][start])+fabs(a[start+1][start+1]))) {
       a[start][start+1] = 0.0;
       start += 1;
       addflops (2);
     }
     if (start >= lim-1) break;
/* Find end of first block */
     end = start+1;
     while (end < lim-1 &&
            fabs(a[end][end+1]) >
            EPS*(fabs(a[end][end])+fabs(a[end+1][end+1]))) {
       end += 1;
       addflops (2);
     }
     if (end < lim-1) a[end][end+1] = 0.0;
     end -= 1;
/* Do shift on upper two by two block */
     u1 = a[start][start];
     u2 = a[start][start+1];
     u3 = a[start+1][start+1];
     l2 = a[end+1][end+1];
       tmp = max (max(fabs(u1),fabs(u2)), fabs(l2));
       u1 /= tmp;
       u2 /= tmp;
       l2 /= tmp;
       r = (u1-l2)*(u1+l2);
       s = u1*u2;
       addflops (7);
       if (r <= 0) { 
	r = 0.0;
	s = 1.0;
       }
       if (fixup (r, s, tmp)) {
         for (k=start; k< min(lim,start+2); k++) {
       	   rotate (a[k][start], a[k][start+1], r, s, tmp);
         }
         if (v != 0) for (k=0; k<cols; k++) {
		  	   rotate (v[start][k], v[start+1][k], r, s, tmp);
         }
         chase (u, a, v, rows, cols, start, end);
       }
   }
   sortsvs (u, a, v, rows, cols);
}


void chase (q, a, v, rows, cols, istart, iend)

/************************************************************************
 *  Purpose:  This subroutine rebidiagonalizes the matrix  	 	*
 *            by chasing unwanted entries down the matrix               *
 *            using orthogonal Givens rotations		     		*
 *            The matrix is factored as q * bidiag * v                  *
 * 									*
 *  Inputs:   a     - (double) array containing matrix to be factored	*
 *            rows  - (integer) number of rows in a		     	*
 *            cols  - (integer) number of columns in a		     	*
 *            istart- (integer) column to start the chasing             *
 *  Outputs:  a     - (double) upper triangular matrix replaces a	*
 *            q     - (double) orthogonal matrix overwrites q		*
 *            v     - (double) orthogonal matrix overwrites v           *
 ************************************************************************
 */									
double **a, **q, **v;
int rows, cols, iend, istart;
{
   int i, ipiv, j, k, lim;
   double r, s, tmp;

   lim = min (rows, cols);
   
   for (i=istart; i<=min(lim,iend); i++) {/*  For each diagonal element */
      ipiv = i;
      j = i+1;
      r = a[ipiv][i];
      s = a[j][i];
      if (fixup (r, s, tmp)) {
         for (k=i; k<min(lim,i+3); k++) { 	/*  Do elimination */
            rotate (a[ipiv][k], a[j][k], r, s, tmp);
         }                           	/*  Closes k loop */ 
         if (q != 0) for (k=0; k<rows; k++) { 	/*  Do multiplication on Q */
            rotate (q[k][ipiv], q[k][j], r, s, tmp);
         }
         a[j][i] = 0.0;  		/*  Kill it */
      }				/*  Closes if */
      ipiv = i+1;
      j = i+2;
      if (j < lim) {
        r = a[i][ipiv];
        s = a[i][j];
        if (fixup (r, s, tmp)) {
          for (k=i; k<(lim,i+3); k++) {    /*  Multiply columns */
            rotate (a[k][ipiv], a[k][j], r, s, tmp);
          }      
          if (v != 0) for (k=0; k<cols; k++) {    /*  Multiply rows of V */
            rotate (v[ipiv][k], v[j][k], r, s, tmp);
          }
          a[i][j] = 0.0;		/*  Kill it */	    
        }				/*  Close if  */
      }
   }					/*  Closes i loop */      
}

void sortsvs (u, s, v, rows, cols)
double **u, **s, **v;
int rows, cols;
{
  double smax, tmp;
  int i,j, lim, imax;
  lim = min (rows, cols);
  for (i=0; i<lim; i++)
    if (s[i][i] < 0.0) {
      if (u != 0) for (j=0; j<rows; j++) u[j][i] = -u[j][i];
      s[i][i] = - s[i][i];
    }

  for (i=0; i<lim; i++) {
    smax = 0.0;
    imax = i;
    for (j=i; j<lim; j++) 
      if (s[j][j] > smax) {	/*  Find maximum singular value */
	smax = s[j][j];
	imax = j;
      }
    if (i != imax) {
      swap (s[imax][imax],s[i][i],tmp);	/* Swap singular values */
      if (u != 0) for (j=0; j<rows; j++) {    
        swap (u[j][imax],u[j][i],tmp);	/* Swap left singular vectors */
      }
      if (v != 0) for (j=0; j<cols; j++) {
        swap (v[imax][j],v[i][j],tmp); /* Swap v rows */
      }
    }
  }
}

m_print (mesg, a, m, n) 
char *mesg;
double **a;
int m, n;
{
  int i,j;
  printf ("\n   %s\n", mesg);
  for (i=0; i<m; i++) {
    for (j=0; j<n; j++) printf ("%10g", a[i][j]);
    printf ("\n");
  }
}
    
void m_ginv (a, b, rows, cols)
double **a, **b;
int rows, cols;
{
	double **u, **v, *s;
	int lim, i, j, k;
	lim = min (rows, cols);
	u = array (rows, rows);
	ident (u, rows, rows);
	v = array (cols, cols);
	ident (v, cols, cols);
	m_svd (u, b, v, rows, cols);
	s = MALLOC (lim, double);
	for (i=0; i<lim; i++) {
		s[i] = (b[i][i] > lim*EPS*b[0][0])? (1.0/b[i][i]) : 0.0;
	}
	addflops (2);
	for (i=0; i<rows; i++) {
		for (j=0; j<cols; j++) {
			a[j][i] = 0.0;
			for (k=0; k<lim; k++) {
				a[j][i] += u[i][k]*s[k]*v[k][j];
			}
			if (lim > 0) addflops (lim);
		}
	}
	FREE (s);
	free_array (u);
	free_array (v);
}

		
	

