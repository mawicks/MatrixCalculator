#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "constant.h"
#include "functions.h"

void dm_mul (d, a, result, rows, cols)
double d, **a, **result;
/*  Multiplies a double times a matrix */
int rows, cols;
{
	int i,j;
	for (i=0; i<rows; i++)
	for (j=0; j<cols; j++)
		result[i][j] = d*a[i][j];
	if (rows>0 && cols>0) addflops (rows*cols);
}

void mm_mul (a, b, c, crows, acols, ccols)
/* Multiplies a matrix times a matrix*/
double **a, **b, **c;
int crows, acols, ccols;
{
	int i,j,k;
	for (i=0; i<crows; i++)
	for (j=0; j<ccols; j++) {
		c[i][j] = 0.0;
		for (k=0; k<acols; k++)
			c[i][j] += a[i][k]*b[k][j];
	}
	if (crows>0 && ccols>0 && acols>0) addflops (2*crows*ccols*acols);
}

void mt_mul (a, b, c, arows, acols, brows)
/* Multiplies a matrix times the transpose of another matrix */
double **a, **b, **c;
int arows, acols, brows;
{
	int i, j, k;
	for (i=0; i<arows; i++)
	for (j=0; j<brows; j++) {
		c[i][j] = 0.0;
		for (k=0; k<acols; k++)
			c[i][j] += a[i][k]*b[j][k]; 
	}
	if (arows>0 && brows>0 && acols>0) addflops (2*arows*brows*acols);
}

void tt_mul (a, b, c, arows, acols, brows)
/* Multiplies the transpose of a matrix times the transpose of a matrix */
double **a, **b, **c;
int arows, acols, brows;
{
	int i, j, k;
	for (i=0; i<acols; i++)
	for (j=0; j<brows; j++) {
		c[i][j] = 0.0;
		for (k=0; k<arows; k++)
			c[i][j] += a[k][i]*b[j][k];
	}
	if (acols>0 && brows>0 && arows>0) addflops (2*acols*brows*arows);
}

void tm_mul (a, b, c, arows, acols, bcols)
/* Multiplies the transpose of a matrix times a matrix */
double **a, **b, **c;
int arows, acols, bcols;
{
	int i, j, k;
	for (i=0; i<acols; i++)
	for (j=0; j<bcols; j++) {
		c[i][j] = 0.0;
		for (k=0; k<arows; k++)
			c[i][j] += a[k][i]*b[k][j]; 
	}
	if (acols>0 && bcols>0 && arows>0) addflops(2*acols*bcols*arows);
}
		
void mm_plus (a, b, c, rows, cols)
/* Adds two matrices */
double **a, **b, **c;
int rows, cols;
{
	int i,j;
	for (i=0; i<rows; i++)
	for (j=0; j<cols; j++)
		c[i][j] = a[i][j] + b[i][j];
	if (rows>0 && cols>0) addflops(rows*cols);
}

void mm_minus (a, b, c, rows, cols)
/* Subtracts two matrices */
double **a, **b, **c;
int rows, cols;
{
	int i,j;
	for (i=0; i<rows; i++)
	for (j=0; j<cols; j++)
		c[i][j] = a[i][j] - b[i][j];
	if (rows>0 && cols>0) addflops(rows*cols);
}

void m_trn (a, b, rows, cols)
/* Transposes a matrix */
double **a, **b;
int rows, cols;
{
	int i,j;
	for (i=0; i<rows; i++)
	for (j=0; j<cols; j++)
		b[j][i] = a[i][j];
}

void m_chs (a, b, rows, cols)
/* Changes the sign of a matrix */
double **a, **b;
int rows, cols;
{
  int i,j;
  for (i=0; i<rows; i++)
    for (j=0; j<cols; j++)
      b[i][j] = - a[i][j];
}

void m_copy (a, b, rows, cols)
/* Copies one matrix to another */
double **a, **b;
int rows, cols;
{
	int i,j;
	for (i=0; i<rows; i++)
	for (j=0; j<cols; j++)
		b[i][j] = a[i][j];
}

double m_rnorm (a, srow, scol, cols)
/* Computes the row norm squared of a row */
double **a;
int srow, scol, cols;
{
	int k;
	double result = 0.0;
	for (k=scol; k<cols; k++) {
		result += a[srow][k] * a[srow][k];
	}
	if (cols > scol) addflops(2*(cols-scol));
	return (result);
}

double m_cnorm (a, srow, scol, rows)
/* Computes the column norm */
double **a;
int srow, scol, rows;
{
	int i;
	double result = 0.0;
	for (i=srow; i<rows; i++) {
		result += a[i][scol] * a[i][scol];
	}
	if (rows > srow) addflops(2*(rows-srow));
	return (result);
}

double m_norm (a, srow, scol, rows, cols)
/* Computes the Frobenius norm of a matrix */
double **a;
int srow, scol, rows, cols;
{
	int i;
	double result = 0.0;
	for (i=srow; i<rows; i++) {
		result += m_rnorm (a, i, scol, cols);
	}
	if (rows > srow) addflops (rows-srow);
	return (result);
}

void rswap (a, cols, r1, r2)
/* Swaps two rows in a matrix */
double **a;
int cols, r1, r2;
{
	int j;
	double tmp;
	for (j=0; j<cols; j++) {
		tmp = a[r1][j];
		a[r1][j] = a[r2][j];
		a[r2][j] = tmp;
	}
}

void cswap (a, rows, c1, c2)
/* Swaps two columns */
double **a;
int rows, c1, c2;
{
	int i;
	double tmp;
	for (i=0; i<rows; i++) {
		tmp = a[i][c1];
		a[i][c1] = a[i][c2];
		a[i][c2] = tmp;
	}
}
	
void ident (a, rows, cols)
/* Creates an identity in a */
double **a;
int rows, cols;
{
	int i,j;
	for (i=0; i<rows; i++) 
	for (j=0; j<cols; j++) a[i][j] = (i==j)? 1.0: 0.0;
}

void m_hes (a, q, size)
/* Makes A upper Hessenberg by a similarity transformation */
double **a, **q;
int size;
{
	int j, k, l, m;
	double cnrm, nrm, qnrm, **v, **tmp1, vmax;
	v = array (1, size);
	tmp1 = array (1, size);
	for (j=0; j<size-1; j++) { /* For each column */
		for (vmax=0.0,k=0; k<size; k++) {
			v[0][k] = k<j+1? 0.0: a[k][j];
			if (fabs(v[0][k]) > vmax) vmax = fabs(v[0][k]);
		}
		if (vmax < SMALL) continue;
		for (k=j+1; k<size; k++) v[0][k] /= vmax;
		if (size>j+1) addflops (size-j-1);
		cnrm = (m_rnorm (v, 0, j+1, size));
		nrm = sqrt (cnrm);
		addsqrts (1);
		qnrm = 2.0 * (cnrm + fabs (v[0][j+1])*nrm);
		v[0][j+1] += (v[0][j+1]>0.0)? nrm: -nrm;
		addflops (4);
		mm_mul (v, a, tmp1, 1, size, size);
		for (l=j+1; l<size; l++) for (m=0; m<size; m++)
			a[l][m] += -2.0 * v[0][l] * tmp1[0][m] / qnrm;
		if (size >j+1) addflops (4*(size-j-1));
		mt_mul (v, q, tmp1, 1, size, size);
		for (l=0; l<size; l++) for (m=j+1; m<size; m++)
			q[l][m] += -2.0 * tmp1[0][l] * v[0][m] / qnrm;
		if (size >0) addflops (4*size);
		mt_mul (v, a, tmp1, 1, size, size);
		for (l=0; l<size; l++) for (m=j+1; m<size; m++)
			a[l][m] += -2.0 * tmp1[0][l] * v[0][m] / qnrm;
		if (size >0) addflops (4*size);
		for (k=j+2; k<size; k++) a[k][j] = 0.0;
	}
	free_array (tmp1);
	free_array (v);
}

void m_shft (a, q, size, idx, eidx)  	/* Do an implicit origin shift */
double **a, **q;			/* using the quadratic given by */
int size, idx, eidx;			/* p and q */
{
	int k, l, m;
	double r, s;
	double cnrm, nrm, qnrm, **v, **tmp1, vmax;
	if (idx+3 > size) {
		printf ("Argument error in m_shift\n");
		exit (1);
	}
	r = a[eidx-1][eidx-1] + a[eidx][eidx];
	s = a[eidx-1][eidx-1]*a[eidx][eidx] -
	    a[eidx-1][eidx]*a[eidx][eidx-1];
	addflops (4);
	v = array (1, size);
	tmp1 = array (1, size);
	for (k=0; k<size; k++) v[0][k] = 0.0;
	v[0][idx] = a[idx][idx]*a[idx][idx]-r*a[idx][idx]+s+
	            a[idx+1][idx]*a[idx][idx+1];
	v[0][idx+1] = a[idx+1][idx]*
		      (a[idx][idx]+a[idx+1][idx+1]-r);
	v[0][idx+2] = a[idx+2][idx+1]*a[idx+1][idx];
	addflops (10);
	for (vmax=0.0,k=idx; k<idx+3; k++) 
		if (fabs(v[0][k]) > vmax) vmax = fabs(v[0][k]);
	if (vmax < SMALL) return;
	for (k=idx; k<idx+3; k++) v[0][k] = v[0][k] / vmax;
	addflops (3);
	cnrm = (m_rnorm (v, 0, idx, idx+3)); /* last parm is 1+last */
	nrm = sqrt (cnrm);
	addsqrts (1);
	qnrm = 2.0 * (cnrm + fabs (v[0][idx])*nrm);
	v[0][idx] += (v[0][idx]>0.0)? nrm: -nrm;
	addflops (4);
	mm_mul (v, a, tmp1, 1, size, size);
	for (l=idx; l<idx+3; l++) for (m=0; m<size; m++)
		a[l][m] += -2.0 * v[0][l] * tmp1[0][m] / qnrm;
	addflops (4*3);
	mt_mul (v, q, tmp1, 1, size, size);
	for (l=0; l<size; l++) for (m=idx; m<idx+3; m++)
		q[l][m] += -2.0 * tmp1[0][l] * v[0][m] / qnrm;
	if (size > 0) addflops (size*4);
	mt_mul (v, a, tmp1, 1, size, size);
	for (l=0; l<size; l++) for (m=idx; m<idx+3; m++)
		a[l][m] += -2.0 * tmp1[0][l] * v[0][m] / qnrm;
	if (size > 0) addflops (size*4);
	free_array (tmp1);
	free_array (v);
}

void m_qrp (a, q, p, rows, cols)
double **a, **q, **p;
int rows, cols;
{
	int j, k, jmax, l, m;
	double maxcnrm, cnrm, nrm, qnrm, **v, **tmp1, **tmp2, vmax;
	v = array (1, rows);
	tmp1 = array (1, cols);
	tmp2 = array (1, rows);
	ident (q, rows, rows);		/* Start with q and p equal to	*/
	ident (p, cols, cols);		/* identies.			*/
	for (j=0; j<cols && j<rows; j++) {	/* For each column, reduce */
		maxcnrm = 0.0;
		jmax = j;
		for (k=j; k<cols; k++) {	/* Find max col norm */
			cnrm = m_cnorm (a, j, k, rows);
			if (cnrm > maxcnrm) {
				maxcnrm = cnrm;
				jmax = k;
			}
		}
		cswap (a, rows, j, jmax);	/* Col with max norm is */
		rswap (p, cols, j, jmax);	/* first now.. 		*/
		for (vmax=0.0,k=0; k<rows; k++) {
			v[0][k] = (k<j)? 0.0: a[k][j];
			if (fabs(v[0][k])>vmax) vmax=fabs(v[0][k]);
		}
		for (k=j; k<rows; k++) v[0][k] /= vmax;
		if (rows > j) addflops (rows-j);
		if (vmax < SMALL) return;
		maxcnrm = m_rnorm (v, 0, j, rows);
		nrm = sqrt (maxcnrm);
		addsqrts (1);
		qnrm = 2.0 * (maxcnrm + fabs(v[0][j])*nrm);
		v[0][j] += (v[0][j] > 0.0)? nrm: -nrm;
		addflops (4);
		mm_mul (v, a, tmp1, 1, rows, cols);
		for (l=j; l<rows; l++) for (m=0; m<cols; m++) 
			a[l][m] += -2.0 * v[0][l] * tmp1[0][m] / qnrm;
		if (rows>j) addflops (4*(rows-j));
		mt_mul (v, q, tmp2, 1, rows, rows);
		for (l=0; l<rows; l++) for (m=j; m<rows; m++)
			q[l][m] += -2.0 * tmp2[0][l] * v[0][m] / qnrm;
		if (rows>0) addflops (4*rows);
		for (k=j+1; k<rows; k++) a[k][j] = 0.0;
	}

	free_array (tmp1);
	free_array (tmp2);
	free_array (v);
}


void m_solve (r, b, rows, bcols)
double **r, **b;
int rows, bcols;
{
	int i, j, k;
	for (i=rows-1; i>=0; i--)
	for (j=0; j<bcols; j++) {		
		for (k=i+1; k<rows; k++)
			b[i][j] -= b[k][j] * r[i][k];
		if (rows > i+1) addflops (2*(rows-i+1));
		b[i][j] /= r[i][i];
		addflops (1);
	}
}

void t_solve (r, b, rows, brows)
double **r, **b;
int rows, brows;
{
	int i, j, k;
	for (i=rows-1; i>=0; i--)
	for (j=0; j<brows; j++) {		
		for (k=i+1; k<rows; k++)
			b[j][i] -= b[j][k] * r[i][k];
		if (rows > i+1) addflops (2*(rows-i+1));
		b[j][i] /= r[i][i];
		addflops (1);
	}
}
