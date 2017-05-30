# MatrixCalculator
Matrix factorizations written in C with an RPN interface

This is some VERY OLD code I wrote in school so I could perform various matrix factorizations and use them interactively in custom algorithms.
It was written before free tools such as Octave, R, etc., were available and even Matlab wasn't available a major universities.
At the time most open-source code for linear algebra existed only in Fortran.
It's probably not very exciting today, but it might be of some interest because the matrix factorizations
(including eigenvalue and singular value decompositions) were written from scratch in C, making it self-contained (no linpack/blas, etc.)

WARNING: The "calculator" interface contains an unsafe gets() call, deprecated today but common when it was written,
which doesn't have any bearing on the matrix factorizations themselves.

The calculator interface is RPN and is modeled after an HP 48G with additional matrix functions (which are listed in "dofunc.c"):

```
$ 1 2 3
$ 5 4 3
$ 2 3 mat
1:  [          1           2           3 ]
    [          5           4           3 ]

$ dup
1:	[          1           2           3 ]
	[          5           4           3 ]

$ qr
1:	[    5.09902     4.31455     3.53009 ]
	[          0     -1.1767    -2.35339 ]

$ *
1:	[          1           2           3 ]
	[          5           4           3 ]

$ svd
1:	[   0.636349    0.575186    0.514024 ]
	[   0.654518   -0.049939   -0.754396 ]
	[  -0.408248    0.816497   -0.408248 ]

$ drop
1:	[    7.77337           0           0 ]
	[          0     1.89068           0 ]

$ drop
1:	[    0.42823    -0.90367 ]
	[    0.90367     0.42823 ]

$ quit
```