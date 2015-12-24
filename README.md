# Krylov Schur algorithm
Krylov Schur eigenvalue algorithm to get eigenvalues/eigenvectors of a matrix.
We do not require this matrix to be formed explicitly; instead, only a function
that returns matrix-vector product is required.

This matlab code is suitable for finding a small subset of eigenvalues/eigenvectors
of a large real matrix. Right now, we do not consider complex matrix case.

#### Reference:
A Krylov-Schur Algorithm for Large Eigenproblems, G. W. Stewart,
SIAM J. M ATRIX A NAL. A PPL ., Vol. 23, No. 3, pp. 601–614
