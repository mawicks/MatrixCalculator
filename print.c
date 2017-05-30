#include <stdio.h>
#include <stdlib.h>

#include "types.h"
#include "constant.h"

void print_item(item)
struct item item;
{
	if (item.datatype == SCALAR) {
		printf ("\t%g\n", item.data.d);
	}
	else if (item.datatype == MATRIX) {
		int i,j;
		for (i=0; i<item.data.mat.rows; i++) {
			printf ("\t[");			
			for (j=0; j<item.data.mat.cols; j++) 
				printf (" %10g ", item.data.mat.m[i][j]);
			printf ("]\n");
		}
	}
	else if (item.datatype == NAME) {
		printf ("\t%s\n", item.data.p);
	}
	else {
		printf ("Internal error in print_item\n");
		exit (1);
	}
}
