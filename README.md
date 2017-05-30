# MatrixCalculator
Matrix factorizations written in C with an RPN interface

This is some very old code I wrote in school so I could perform various matrix factorizations.
It was written before free tools such as Octave, R, etc., were available.  At the time most code for linear algebra exited only in Fortran.
It's probably not very exciting today, but it might be of some interest because the matrix factorizations
(including eigenvalue and singular value decompositions) were written from scratch in C, making it self-contained (no linpack/blas, etc.)

WARNING: The "calculator" interface contains a gets() call, eprecated today but common when it was written,
which doesn't have any bearing on the matrix factorizations themselves.






