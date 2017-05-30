#include "types.h"
#include "functions.h"
#include "constant.h"
#include "mem.h"
#include <stdio.h>
#define NUMFUNCTS	49
#define MAXNAMES	100
typedef void (*ptof) ();

ptof fn_lookup ();
struct item nm_lookup();

static char *fn_table[] = {	"*",		"/",
				"+",		"-",
				"alg1",         "alg2",
				"alg3",         "alg4",
				"drop",		"dr", 
				"dr", 
				"dup",
				"cont1",	"cont2",
				"ev",
				"fnrm",
				"ginv",	"gi",
				"hes",
				"idn",		"ident",
				"inv",		"inverse",
				"mat",		"matrix",
				"over",
				"print",	"pr",
				"pr", 		
				"svd",
				"swap",		"trn",
				"qu",		"quit",
				"qr",
				"qrp", 		"store", 
				"sto",
				"cos",          "sin",
				"tan",
				"exp",          "atan",
				"ln",           "log",
				"flops",        "noflops",
				"sqrt",         "chs"
			  };

static void (*f_table[])() = {	times,		divide,
				plus,		minus,
				alg1,           alg2,
				alg3,           alg4,
				drop,		drop,
				drop,
				dupl,
				cont1,		cont2,
				ev,
				fnrm,
				ginv,		ginv,
				hes,
				idn,		idn,
				inverse,	inverse,
				get_mat,	get_mat,
				over,
				print,		print,
				print,
				svd, 
				swap,		trn,
				quit,		quit,
				qr,
				qrp, 		store,
				store,
				ucos,           usin,
				utan,
				uexp,           uatan,
				uln,            ulog,
				flops,          noflops,
				usqrt,          uchs

			  };
static int 	numnames = 0;
static char	*names[MAXNAMES];
struct item saved_items[MAXNAMES];

void dofunct(name)
char *name;
{

	void (*func)();
	struct item lookup, new;
	if ((func = fn_lookup(name)) != NULL)
		(*func)();
	else if (*name != '\'' && 
	        (lookup=nm_lookup(name), lookup.datatype != EMPTY ))
		push (lookup);
	else {
		if (*name == '\'') name++;
		new.datatype = NAME;
		new.data.p = MALLOC (strlen(name)+1, char);
		strcpy (new.data.p, name);
		push (new);
	}
}


ptof fn_lookup (name)
char *name;
{
	int i;
	for (i=0; i<NUMFUNCTS && strcmp(name, fn_table[i]) != 0; i++);
	if (i <  NUMFUNCTS) {
		return (f_table[i]);
	}
	else {
		return (NULL);
	}
}

struct item nm_lookup (name)
char *name;
{
	int i;
	struct item empty;
	for (i=0; i<numnames && strcmp(name, names[i])!=0; i++);
	if (i < numnames) {
		return (item_copy(saved_items[i]));
	}
	else {
		empty.datatype = EMPTY;
		return (empty);
	}
}

void store_name (name, item)
char *name;
struct item item;
{
	int i;
	for (i=0; i<numnames && strcmp (name, names[i]) !=0; i++);
	if (i< numnames) {
		del_item (saved_items[i]);
		free (names[i]);
		saved_items[i] = item;
		names[i] = name;
	}
	else if (numnames < MAXNAMES) {
		saved_items[numnames] = item;
		names[numnames++] = name;
	}
	else printf ("Too many names...\n");
}
