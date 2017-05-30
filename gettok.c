#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include "mem.h"
#include "types.h"
#include "constant.h"

char *prmpt = "$ ";
double getnum();
char *getident();
char ib[INBUFSIZE] = {'\0'};
char *p = ib;

struct token gettok()
{
	struct token result;
	while (p != NULL && ((*p == ' ') || (*p == '\t') || (*p == '\0'))) {
		if (*p == '\0') { 
			printf (prmpt);
			p = gets(ib);
		}
		else p++;
	}
	if (p != NULL) {
		while ((*p == ' ') || (*p == '\t')) p++;
		if (isdigit (*p) || (*p == '.') ||
		   ((*p) == '-' && isdigit (*(p+1))) ||
                   ((*p) == '-' && *(p+1) == '.')) {
			result.token = SCALAR;
			result.value.d = getnum();
		}
		else if (isalpha (*p) || *p == '\'') {
			result.token = FUNCTION;
			result.value.p = getident();
		}
		else {
			result.token = FUNCTION;
			result.value.p = MALLOC (2, char);
			result.value.p[0] = *p++;
			result.value.p[1] = '\0';
		}
	}
	else {
		result.token = EOF;
	}
	return (result);
}


char *getident ()
{
	char *result, *start = p, *index;
	while (isalnum (*++p));		/* Find end */
	result = MALLOC (p - start + 1, char);	/* Get enough space */
	for (index = result; start<p; index++, start++)
		*index = *start;
	*index = '\0';
	return (result);
}

double getnum()
{
	double result = 0, power = 1.0;
	int digit, sign=1;
	if (*p == '-') {
		p++;
		sign = -1;
	}
	while (isdigit (*p)) {
		digit = *(p++) - '0';	/* Force this first */
		result = 10.0 * result + digit;
	}
	if (*p == '.') {
		p++;
		while (isdigit (*p)) {
			power *= 10.0;
			digit = *(p++) - '0';
			result = result + digit / power;
		}
	}
	result *= sign;
	sign = 1;
	if (*p == 'e' || *p == 'E') {
		p++;
		if (*p == '-') {
			p++;
			sign = -1;
		}
		else if (*p == '+') p++;
		power = 0.0;
		while (isdigit (*p)) {
			digit = *(p++) - '0';
			power = 10.0 * power + digit;
		}
		power *= sign;
		result = result * pow(10.0, power);
	}
	return (result);
}
