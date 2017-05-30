#ifndef mem_h
#include <stdlib.h>

#define mem_h
#define NIL 0
char *Calloc (), *Malloc ();
char *xcalloc (), *xmalloc ();
void clear ();

#define CALLOC(n, type) \
	((type *) xcalloc ((n), sizeof (type)))

#define MALLOC(n, type) \
	((type *) xmalloc ((n)*sizeof (type)))

#define CLEAR(pointer, n, type) \
	( clear ( (pointer), (n)*sizeof (type)))

#define FREE(p) (free(p))
#endif
