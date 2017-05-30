#define TYPES_H

struct matrix { 
	double **m; 
	int rows, cols;
	};

union data {
	struct matrix	mat;
	double		d;
	char		*p;
	};

struct item {
	int	datatype;	/* Identifies matrix or scalar */
	union data data;
	};

union token_value {
	double	d;
	char	*p;
	};

struct token {
	int	token;
	union token_value	value;
};
