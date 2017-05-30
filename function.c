#include "constant.h"
#include "types.h"
#include "functions.h"
#include "math.h"

void print ()
{
	print_stack();
}
void flops ()
{
        onflops();
}
void noflops ()
{
        offflops();
}

void drop ()
{
	struct item item;
	item = pop ();
	del_item (item);
}

void get_mat()
{
	struct item row_item, col_item, entry, matrix;
	int error = 0;
	int rows, cols, i, j;
	mark_stack();
	col_item = pop();
	row_item = pop();
	if (row_item.datatype != SCALAR || 
	    col_item.datatype != SCALAR) {
		printf ("Argument error.\n");
		reset_stack();
	}
	else {
		rows = row_item.data.d;
		cols = col_item.data.d;
		matrix.datatype = MATRIX;
		matrix.data.mat.m = array (rows, cols);
		matrix.data.mat.rows = rows;
		matrix.data.mat.cols = cols;

		for (i=rows-1; !error && i>=0; i--)
		for (j=cols-1; !error && j>=0; j--) {
			if (entry = pop(), entry.datatype == SCALAR)
				matrix.data.mat.m[i][j] = entry.data.d;
			else error = 1;
		}
		if (error) {
			printf ("Argument error.\n");
			reset_stack();
			free_array(matrix.data.mat.m);
		}
		else {
			push (matrix);
		}
	}
}

void times()
{
	struct item a, b, c;
	mark_stack();
	b = pop();
	a = pop();
	if ((a.datatype == SCALAR) && (b.datatype == SCALAR)) {
		c.datatype = SCALAR;
		c.data.d = a.data.d * b.data.d;
		addflops (1);
		push (c);
	}
	else if ((a.datatype == SCALAR) &&
	         (b.datatype == MATRIX)) {
		c = b;
		c.data.mat.m = array(c.data.mat.rows, c.data.mat.cols);
		dm_mul (a.data.d, b.data.mat.m, c.data.mat.m,
		        c.data.mat.rows, c.data.mat.cols);
		free_array (b.data.mat.m);
		push (c);
	}
	else if ((a.datatype == MATRIX) &&
	         (b.datatype == SCALAR)) {
		c = a;
		c.data.mat.m = array(c.data.mat.rows, c.data.mat.cols);
		dm_mul (b.data.d, a.data.mat.m, c.data.mat.m,
		        c.data.mat.rows, c.data.mat.cols);
		free_array (a.data.mat.m);
		push (c);
	}
	else if ((a.datatype == MATRIX) &&
	         (b.datatype == MATRIX) &&
		 (a.data.mat.cols == b.data.mat.rows)) {
		c.datatype = MATRIX;
		c.data.mat.rows = a.data.mat.rows;
		c.data.mat.cols = b.data.mat.cols;
		c.data.mat.m = array(c.data.mat.rows, c.data.mat.cols);
		mm_mul (a.data.mat.m, b.data.mat.m, c.data.mat.m,
		        c.data.mat.rows, a.data.mat.cols, 
		        c.data.mat.cols);
		free_array (a.data.mat.m);
		free_array (b.data.mat.m);
		push (c);
	}
	else {
		reset_stack ();		/* Something must be wrong */
		printf ("Argument error\n");
	}
}

void divide ()
{
	struct item a, b, c;
	mark_stack ();
	b = pop();
	a = pop();
	if ((a.datatype == SCALAR) &&
	    (b.datatype == SCALAR)) {
		c.datatype = SCALAR;
		c.data.d = a.data.d / b.data.d;
		addflops(1);
		push (c);
	    }
	else if ((a.datatype == MATRIX) &&
	         (b.datatype == SCALAR)) {
		c = a;
		c.data.mat.m = array(c.data.mat.rows, c.data.mat.cols);
		dm_mul (1.0/b.data.d, a.data.mat.m, c.data.mat.m,
		        c.data.mat.rows, c.data.mat.cols);
		free_array (a.data.mat.m);
		push (c);
	}
	else if ((a.datatype == MATRIX) &&
		 (b.datatype == MATRIX) &&
		 (b.data.mat.rows == b.data.mat.cols) &&
		 (a.data.mat.rows == b.data.mat.rows)) {
		double **q, **p;
		q = array (b.data.mat.rows, b.data.mat.rows);
		p = array (b.data.mat.cols, b.data.mat.cols);
		c.data.mat.m = array (a.data.mat.rows, a.data.mat.cols);
		m_qrp (b.data.mat.m, q, p, b.data.mat.rows, b.data.mat.cols);
		tm_mul (q, a.data.mat.m, c.data.mat.m, 
		        b.data.mat.rows, b.data.mat.rows, a.data.mat.cols);
		m_solve (b.data.mat.m, c.data.mat.m, b.data.mat.rows, 
		         a.data.mat.cols);
		tm_mul (p, c.data.mat.m, a.data.mat.m,
			b.data.mat.cols, b.data.mat.cols, a.data.mat.cols);
		push (a);
		free_array (q);
		free_array (p);
		free_array (c.data.mat.m);
		free_array (b.data.mat.m);
	}		
	else {
		reset_stack();
		printf ("Argument error.\n");
	}
}

void inverse()
{
	struct item a, b;
	mark_stack();
	a = pop();
	if (a.datatype == SCALAR) {
		a.data.d = 1.0 / a.data.d;
		addflops(1);
		push (a);
	}
	else if (a.datatype == MATRIX &&
		 a.data.mat.rows == a.data.mat.cols) {
		double **q, **p;
		q = array (a.data.mat.rows, a.data.mat.rows);
		p = array (a.data.mat.cols, a.data.mat.cols);
		m_qrp (a.data.mat.m, q, p, a.data.mat.rows, a.data.mat.cols);
		t_solve (a.data.mat.m, q,
		         a.data.mat.rows, a.data.mat.rows);
		b = a;
		b.data.mat.m = array (a.data.mat.rows, a.data.mat.cols);
		tt_mul (p, q, b.data.mat.m, a.data.mat.rows, a.data.mat.cols, 
		        a.data.mat.cols);
		free_array (q);
		free_array (p);
		free_array (a.data.mat.m);
		push (b);
	}
	else {
		reset_stack ();
		printf ("Argument error.\n");
	}
}
void ginv()
{
	struct item a, b;
	mark_stack();
	a = pop();
	if (a.datatype == SCALAR) {
		a.data.d = 1.0 / a.data.d;
		addflops(1);
		push (a);
	}
	else if (a.datatype == MATRIX) {
		b.datatype = MATRIX;
		b.data.mat.m = array (a.data.mat.cols, a.data.mat.rows);
		b.data.mat.rows = a.data.mat.cols;
		b.data.mat.cols = a.data.mat.rows;
		m_ginv (b.data.mat.m, a.data.mat.m, a.data.mat.rows, a.data.mat.cols);
		free_array (a.data.mat.m);
		push (b);
	}
	else {
		reset_stack ();
		printf ("Argument error.\n");
	}
}
void uchs()
{
	struct item a, b;
	mark_stack();
	a = pop();
	if (a.datatype == SCALAR) {
		a.data.d = - a.data.d;
		addflops(1);
		push (a);
	}
	else if (a.datatype == MATRIX) {
		m_chs (a.data.mat.m, a.data.mat.m, a.data.mat.rows, a.data.mat.cols);
		push (a);
	}
	else {
		reset_stack ();
		printf ("Argument error.\n");
	}
}
		
void plus()
{
	struct item a, b, c;
	mark_stack();
	b = pop();
	a = pop();
	if ((a.datatype == SCALAR) && (b.datatype == SCALAR)) {
		c.datatype = SCALAR;
		c.data.d = a.data.d + b.data.d;
		addflops(1);
		push (c);
	}
	else if ((a.datatype == MATRIX) && (b.datatype == MATRIX) &&
	         (a.data.mat.rows == b.data.mat.rows) && 
	         (a.data.mat.cols == b.data.mat.cols)) {
		c = a;
		c.data.mat.m = array(c.data.mat.rows, c.data.mat.cols);
		mm_plus (a.data.mat.m, b.data.mat.m, c.data.mat.m,
		         a.data.mat.rows, a.data.mat.cols);
		free_array (a.data.mat.m);
		free_array (b.data.mat.m);
		push (c);
	}
	else {
		reset_stack();
		printf ("Argument error.\n");
	}
}

		
void minus()
{
	struct item a, b, c;
	mark_stack();
	b = pop();
	a = pop();
	if ((a.datatype == SCALAR) && (b.datatype == SCALAR)) {
		c.datatype = SCALAR;
		c.data.d = a.data.d - b.data.d;
		addflops(1);
		push (c);
	}
	else if ((a.datatype == MATRIX) && (b.datatype == MATRIX) &&
	         (a.data.mat.rows == b.data.mat.rows) && 
	         (a.data.mat.cols == b.data.mat.cols)) {
		c = a;
		c.data.mat.m = array(c.data.mat.rows, c.data.mat.cols);
		mm_minus (a.data.mat.m, b.data.mat.m, c.data.mat.m,
		          a.data.mat.rows, a.data.mat.cols);
		free_array (a.data.mat.m);
		free_array (b.data.mat.m);
		push (c);
	}
	else {
		reset_stack();
		printf ("Argument error.\n");
	}
}


void trn ()
{
	struct item a, b;
	mark_stack();
	a = pop();
	if (a.datatype == SCALAR) push (a);
	else if (a.datatype == MATRIX) {
		b.datatype = MATRIX;
		b.data.mat.rows = a.data.mat.cols;
		b.data.mat.cols = a.data.mat.rows;
		b.data.mat.m = array(b.data.mat.rows, b.data.mat.cols);
		m_trn (a.data.mat.m, b.data.mat.m, 
		       a.data.mat.rows, a.data.mat.cols);
		free_array (a.data.mat.m);
		push (b);
	}
	else {
		reset_stack();
		printf ("Argument error.\n");
	}
}

void usin ()
{
  struct item a;
  mark_stack ();
  a = pop ();
  if (a.datatype == SCALAR) {
    a.data.d = sin (a.data.d);
    push (a);
  }
  else {
    reset_stack();
    printf ("Argument error.\n");
  }
}
void ucos ()
{
  struct item a;
  mark_stack ();
  a = pop ();
  if (a.datatype == SCALAR) {
    a.data.d = cos (a.data.d);
    push (a);
  }
  else {
    reset_stack();
    printf ("Argument error.\n");
  }
}
void usqrt ()
{
  struct item a;
  mark_stack ();
  a = pop ();
  if (a.datatype == SCALAR) {
    a.data.d = sqrt (a.data.d);
    push (a);
  }
  else {
    reset_stack();
    printf ("Argument error.\n");
  }
}
void utan ()
{
  struct item a;
  mark_stack ();
  a = pop ();
  if (a.datatype == SCALAR) {
    a.data.d = tan (a.data.d);
    push (a);
  }
  else {
    reset_stack();
    printf ("Argument error.\n");
  }
}
void uexp ()
{
  struct item a;
  mark_stack ();
  a = pop ();
  if (a.datatype == SCALAR) {
    a.data.d = exp (a.data.d);
    push (a);
  }
  else {
    reset_stack();
    printf ("Argument error.\n");
  }
}

void uatan ()
{
  struct item a;
  mark_stack ();
  a = pop ();
  if (a.datatype == SCALAR) {
    a.data.d = atan (a.data.d);
    push (a);
  }
  else {
    reset_stack();
    printf ("Argument error.\n");
  }
}
void uln ()
{
  struct item a;
  mark_stack ();
  a = pop ();
  if (a.datatype == SCALAR) {
    a.data.d = log (a.data.d);
    push (a);
  }
  else {
    reset_stack();
    printf ("Argument error.\n");
  }
}
void ulog ()
{
  struct item a;
  mark_stack ();
  a = pop ();
  if (a.datatype == SCALAR) {
    a.data.d = log10 (a.data.d);
    push (a);
  }
  else {
    reset_stack();
    printf ("Argument error.\n");
  }
}


void dupl ()
{
	struct item a, b;
	mark_stack();
	a = pop();
	if (a.datatype == SCALAR) {
		push (a);
		push (a);
	}
	else if (a.datatype == MATRIX) {
		push (a);
		b = a;
		b.data.mat.m = array(b.data.mat.rows, b.data.mat.cols);
		m_copy (a.data.mat.m, b.data.mat.m,
		        b.data.mat.rows, b.data.mat.cols);
		push (b);
	}
	else {
		reset_stack();
		printf ("Argument error.\n");
	}
}

void over ()
{
	struct item a, b, c;
	mark_stack();
	b = pop();
	a = pop();
	if (a.datatype == SCALAR) {
		push (a);
		push (b);
		push (a);
	}
	else if (a.datatype == MATRIX) {
		push (a);
		push (b);
		c = a;
		c.data.mat.m = array(a.data.mat.rows, a.data.mat.cols);
		m_copy (a.data.mat.m, c.data.mat.m,
		        c.data.mat.rows, c.data.mat.cols);
		push (c);
	}
	else {
		reset_stack();
		printf ("Argument error.\n");
	}
}

void swap ()
{
	struct item a, b;
	mark_stack ();

	b = pop ();
	a = pop ();
	if ((a.datatype != EMPTY) && (b.datatype != EMPTY)) {
		push (b);
		push (a);
	}
	else {
		reset_stack ();
		printf ("Argument error.\n");
	}
}

void quit ()
{
	exit (0);
}

void hes ()
{
	struct item q, r;
	mark_stack ();
	r = pop ();
	q = pop ();
	if (r.datatype == MATRIX &&
	    r.data.mat.rows == r.data.mat.cols &&
	    q.datatype == MATRIX && 
	    q.data.mat.rows == q.data.mat.cols &&
	    q.data.mat.rows == r.data.mat.rows) {
		m_hes (r.data.mat.m, q.data.mat.m,
		       r.data.mat.rows);
		push (q);
		push (r);
	}
	else {
		printf ("Argument error.\n");
		reset_stack ();
	}
}
		
void ev ()
{
	struct item q, r;
	mark_stack ();
	r = pop();

	if (r.datatype == MATRIX &&
	    r.data.mat.rows == r.data.mat.cols) {
		q = r;
		q.data.mat.m=array (q.data.mat.rows, q.data.mat.cols);
		eigv (r.data.mat.m, q.data.mat.m, r.data.mat.rows);
		push (q);
		push (r);
	}
	else {
		printf ("Argument error.\n");
		reset_stack ();
	}
}
		
void qr ()
{
	struct item q, r;
	mark_stack ();
	r = pop();
	if (r.datatype == MATRIX) {
		q.datatype = MATRIX;
		q.data.mat.rows = r.data.mat.rows;
		q.data.mat.cols = r.data.mat.rows;
		q.data.mat.m = array (r.data.mat.rows, r.data.mat.rows);
		m_qr (r.data.mat.m, q.data.mat.m,
		       r.data.mat.rows, r.data.mat.cols);
		push (q);
		push (r);
	}
	else {
		printf ("Argument error.\n");
		reset_stack ();
	}
}
		
void qrp ()
{
	struct item q, r, p;
	mark_stack ();
	r = pop();
	if (r.datatype == MATRIX) {
		q.datatype = MATRIX;
		q.data.mat.rows = r.data.mat.rows;
		q.data.mat.cols = r.data.mat.rows;
		p.datatype = MATRIX;
		p.data.mat.rows = r.data.mat.cols;
		p.data.mat.cols = r.data.mat.cols;
		q.data.mat.m = array (r.data.mat.rows, r.data.mat.rows);
		p.data.mat.m = array (r.data.mat.cols, r.data.mat.cols);
		m_qrp (r.data.mat.m, q.data.mat.m, p.data.mat.m,
		     r.data.mat.rows, r.data.mat.cols);
		push (q);
		push (r);
		push (p);
	}
	else {
		printf ("Argument error.\n");
		reset_stack ();
	}
}
		
void idn ()
{
	struct item a, b;
	mark_stack ();
	a = pop();
	if (a.datatype == SCALAR) {
		b.datatype = MATRIX;
		b.data.mat.rows = a.data.d;
		b.data.mat.cols = a.data.d;
		b.data.mat.m = array (b.data.mat.rows, b.data.mat.cols);
		ident (b.data.mat.m, b.data.mat.rows, b.data.mat.cols);
		push (b);
	}
	else {
		reset_stack ();
		printf ("Argument error.\n");
	}
}

void cont1 ()
{
	struct item a, b, q, norm;
	mark_stack();
	a = pop();
	b = pop();
	q = pop();
	if (a.datatype == MATRIX &&
	    b.datatype == MATRIX &&
	    q.datatype == MATRIX &&
	    a.data.mat.rows == a.data.mat.cols &&
	    b.data.mat.rows == a.data.mat.rows &&
	    q.data.mat.rows == q.data.mat.cols &&
	    q.data.mat.rows == a.data.mat.rows) {
	        norm.datatype = SCALAR;
		control1 (a.data.mat.m, b.data.mat.m, q.data.mat.m,
		         a.data.mat.rows, b.data.mat.cols, &norm.data.d);
		push (q);
		push (b);
		push (a);
		push (norm);
	    }
	else {
		printf ("Argument error.\n");
		reset_stack();
	}
}
void alg1 ()
{
	struct item a, b, q, norm;
	mark_stack();
	a = pop();
	b = pop();
	q = pop();
	if (a.datatype == MATRIX &&
	    b.datatype == MATRIX &&
	    q.datatype == MATRIX &&
	    a.data.mat.rows == a.data.mat.cols &&
	    b.data.mat.rows == a.data.mat.rows &&
	    q.data.mat.rows == q.data.mat.cols &&
	    q.data.mat.rows == a.data.mat.rows) {
	        norm.datatype = SCALAR;
		algo1 (q.data.mat.m, b.data.mat.m, a.data.mat.m,
 		       a.data.mat.rows, b.data.mat.cols, &norm.data.d);
		push (q);
		push (b);
		push (a);
		push (norm);
	    }
	else {
		printf ("Argument error.\n");
		reset_stack();
	}
}
void alg2 ()
{
	struct item a, b, q, norm;
	mark_stack();
	a = pop();
	b = pop();
	q = pop();
	if (a.datatype == MATRIX &&
	    b.datatype == MATRIX &&
	    q.datatype == MATRIX &&
	    a.data.mat.rows == a.data.mat.cols &&
	    b.data.mat.rows == a.data.mat.rows &&
	    q.data.mat.rows == q.data.mat.cols &&
	    q.data.mat.rows == a.data.mat.rows) {
	        norm.datatype = SCALAR;
		algo2 (q.data.mat.m, b.data.mat.m, a.data.mat.m,
 		       a.data.mat.rows, b.data.mat.cols, &norm.data.d);
		push (q);
		push (b);
		push (a);
		push (norm);
	    }
	else {
		printf ("Argument error.\n");
		reset_stack();
	}
}
void alg3 ()
{}
void alg4 ()
{}

void cont2 ()
{
	struct item a, b, q, norm;
	mark_stack();
	a = pop();
	b = pop();
	q = pop();
	if (a.datatype == MATRIX &&
	    b.datatype == MATRIX &&
	    q.datatype == MATRIX &&
	    a.data.mat.rows == a.data.mat.cols &&
	    b.data.mat.rows == a.data.mat.rows &&
	    q.data.mat.rows == q.data.mat.cols &&
	    q.data.mat.rows == a.data.mat.rows) {
	        norm.datatype = SCALAR;
		control2 (a.data.mat.m, b.data.mat.m, q.data.mat.m,
		         a.data.mat.rows, b.data.mat.cols,&norm.data.d);
		push (q);
		push (b);
		push (a);
		push (norm);
	    }
	else {
		printf ("Argument error.\n");
		reset_stack();
	}
}


void fnrm ()
{
	struct item a, b;
	mark_stack ();
	a = pop();
	if (a.datatype == MATRIX) {
		b.datatype = SCALAR;
		b.data.d = sqrt (m_norm(a.data.mat.m, 0, 0, 
		                 a.data.mat.rows, a.data.mat.cols));
		free_array (a.data.mat.m);
		push (b);
	}
	else {
		reset_stack ();
		printf ("Argument error.\n");
	}
}

void store ()
{
	struct item a, nm;
	mark_stack ();
	nm = pop();
	a = pop();
	if (nm.datatype == NAME && a.datatype != EMPTY) {
		store_name (nm.data.p, a);
	}
	else  {
		printf ("Argument error.\n");
		reset_stack();
	}
}

void svd ()
{
	struct item u, s, v;
	mark_stack ();
	s = pop ();
	if (s.datatype == MATRIX) {
		u.datatype = MATRIX;
		u.data.mat.rows = s.data.mat.rows;
		u.data.mat.cols = s.data.mat.rows;
		u.data.mat.m = array (u.data.mat.rows, u.data.mat.cols);
		v.datatype = MATRIX;
		v.data.mat.rows = s.data.mat.cols;
		v.data.mat.cols = s.data.mat.cols;
		v.data.mat.m = array (v.data.mat.rows, v.data.mat.cols);
		m_svd (u.data.mat.m, s.data.mat.m,
		       v.data.mat.m, s.data.mat.rows, s.data.mat.cols);
		push (u);
		push (s);
		push (v);
	}
	else {
		printf ("Argument error.\n");
		reset_stack();
	}
}
	
