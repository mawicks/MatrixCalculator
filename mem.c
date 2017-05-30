/* Error checking memory management routines. */
#include <stdio.h>
#include <stdlib.h>

#include "mem.h"

static char *message = "Error allocating memory!\n";

char *xcalloc (n, size) 
int n, size;
{
	char *p;
	if ((p = calloc (n, size)) == NIL) {
		printf (message);
		exit (1);
	}
	return (p);
}

char *xmalloc (size)
int size;
{
	char *p;
	if ((p = malloc (size)) == NIL) {
		printf (message);
		exit (1);
	}
	return (p);
}

void clear (area, size) 
register char *area;
int size;
{
	for (; size > 0; size--) *(area++) = 0;
}		
