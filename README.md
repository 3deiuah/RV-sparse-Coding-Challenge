# RV-sparse-Coding-Challenge
Coding challenge for RV-sparse under the LFX mentorship program
The function assumes all memory has been preallocated
There are 2 components to the ```sparse_multiply``` function, a dense to sparse matrix parser, and a sparse matrix multiply
# Dense Matrix to CSR Sparse Matrix
The Dense Matrix to CSR Sparse Matrix assumes a row major, row x cols matrix A such that
$$A_{ij}=$$```A[i*cols+j]```
Given the dimensions of the matrix by the ```rows``` and ```cols``` variables, $$i*\text{cols}$$ will always index the first element of the $$i$$-th row(i.e. $$A_{11},A_{21},A_{31}...$$). Constraining $$0 \leq i \leq$$rows ensures that A will not index out of bounds (i.e. $$A_{rows+1,1}$$)
Then constraining $$0 \leq j \leq$$cols ensures A will not index out of bounds column wise. 
This does not account for the case in which the arguments are poorly defined
# CSR Sparse Matrix multiplication
The sparse matrix multiplication is done iteratively over the rows of the matrix, then iteratively over the nonzero values contained in the row.
# Running the program
The project can be run by entering:
``` gcc -lm -o run challenge.c
./run
```