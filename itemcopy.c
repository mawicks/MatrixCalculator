#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "functions.h"
#include "mem.h"
#include "types.h"
#include "constant.h"

struct item item_copy (item)
struct item item;
{
	struct item result;
	int i,j;
	if (item.datatype == SCALAR) {
		result = item;
	}
	else if (item.datatype == MATRIX) {
		result.data.mat.rows = item.data.mat.rows;
		result.data.mat.cols = item.data.mat.cols;
		result.data.mat.m = array (result.data.mat.rows, 
				           result.data.mat.cols);
		result.datatype = MATRIX;
		
		for (i=0; i<result.data.mat.rows; i++)
		for (j=0; j<result.data.mat.cols; j++)
			result.data.mat.m[i][j] = item.data.mat.m[i][j];
	}
	else if (item.datatype == NAME) {
		result = item;
		result.data.p = MALLOC (strlen(item.data.p)+1, char);
		strcpy (result.data.p, item.data.p);
	}
	else {
		printf ("Internal error in item_copy\n");
		exit (1);
	}
	return (result);
}


void del_item (item)
struct item item;
{
	if (item.datatype == MATRIX) {
		free_array (item.data.mat.m);
	}
	else if (item.datatype == NAME) {
		free(item.data.p);
		}
}
