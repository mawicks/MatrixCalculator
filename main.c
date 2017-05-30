#include <stdio.h>
#include "constant.h"
#include "types.h"
#include "functions.h"

/* Main routine for matrix calculater */

main ()
{
	struct token token;
	struct item item;
	while (token = gettok(), token.token != EOF) {
		if (token.token == SCALAR) {
			item.datatype = SCALAR;
			item.data.d = token.value.d;
			push (item);
		}
		else if (token.token == FUNCTION) {
			dofunct(token.value.p);
			free (token.value.p);
			print_top ();
			prtflops ();
		}
		else if (token.token == ERROR) {
			printf ("Input error - ignored\n");
		}
		else {			/* This should not happen */
			printf ("Internal error in main!\n");
			exit (1);
		}
	}
}
