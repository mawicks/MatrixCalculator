/*
 * These routines manage a stack of scalars or matrices for a simple stack
 * oriented matrix/scalar calculater
 */
#include <stdio.h>

#include "mem.h"
#include "types.h"
#include "constant.h"
#include "functions.h"

static struct item stack[STACKSIZE];
static int sp = 0;

void push (item)
struct item item;
{
	if (sp == STACKSIZE) {
		printf ("Stack limit exceeded\n");
		exit (1);
	}
	else if ((item.datatype == SCALAR) || 
	         (item.datatype == MATRIX) ||
		 (item.datatype == NAME)) {
		stack[sp++] = item;
	}
	else {				/* Can't happen */
		printf ("Internal inconsistancy in push\n");
		exit (1);
	}
}

struct item pop ()
{
	struct item value;
	if (sp == 0) {
		printf ("Stack is empty!\n");
		value.datatype = EMPTY;
	}
	else {
		value = stack[--sp];
	}
	return (value);
}


void print_stack()
{
	int i;
	for (i=0; i<sp-1; i++) {
		printf ("%d: ", sp-i);
		print_item (stack[i]);
		printf ("\n");
	}
}


void print_top()
{
	if (sp != 0) {
	printf ("%d: ", 1);
	print_item (stack[sp-1]);
	printf ("\n");
	}
}

static int mark;
void mark_stack()
{
	mark = sp;
}

void reset_stack()
{
	sp = mark;
}
