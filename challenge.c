#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// =========================================================
// FUNCTION PROTOTYPE
// =========================================================
void sparse__multiply(
    int rows,
    int cols,
    const double* A,
    const double* x,
    int* out_nnz,
    double* values,
    int* col_indices,
    int* row_ptrs,
    double* y
);
// Let A be a row x cols matrix 
// Assume row major(i.e. A[n_1] is on the same row as A[n_2] IFF n_1 and n_2 < m*cols for the same integer m)
// Let col_indicies be the column indicies of the non-zero values of A
// Let row_ptrs be the row index range of the non-zero values of A
// Let values be the non-zero values of A
// Let out_nnz be an integer of the number of non-zero numbers in A
// Let x be a cols x 1 matrix
// Find y = A*x(where y is a rows x 1 matrix)



void sparse_multiply(
    int rows, int cols, const double* A, const double* x,
    int* out_nnz, double* values, int* col_indices, int* row_ptrs,
    double* y
) {
    //Scans a row-major matrix A and identifies its non-zero elements.
    //Extracts them into Compressed Sparse Row (CSR) format using caller-provided buffers.


    //Initialize temporary number of non-zero numbers as 0
    int current_nnz= 0;

    //Iterate over rows
    for (int i =0; i<rows;i++) {

        //Assign row_ptrs[i] to the number of non-zero numbers in the previous row
        row_ptrs[i]=current_nnz;

        //Get value of A_{ij} 
        //ASSUME ROW-MAJOR
        //row i of A will contain cols number of values, so A[i*cols] will retrive the first element of the i-th row(i.e. A_{11}, A_{11}, A_{31}...)
        int current_row = i*cols;

        //Iterate over columns
        //j = current column index
        for (int j = 0; j < cols;j++) {

            //To retrieve the j-th element of the row, add j(<cols) to the current_row value
            double val = A[current_row+j];

            //If val !=0
            if (val) {
                //Append current value to list of non-zero values
                values[current_nnz]=val;
                //Append current column index to column indicies
                col_indices[current_nnz] = j;
                //Iterate number of non-zero numbers
                current_nnz++;
            }
        }
    }
    //Assign the last value of row_ptrs to the number of non-zero values
    row_ptrs[rows]=current_nnz;
    *out_nnz = current_nnz;
    //Computes the matrix-vector product y = A * x using the extracted CSR data.
    //Writes the result directly into a caller-provided output buffer.

    for (int i = 0;i < rows; i++) {
        //Initialize default output to zero
        y[i]=0.0;
        //Get the range of the non-zero values from row_ptrs
        int row_start = row_ptrs[i];
        int row_end = row_ptrs[i+1];
        //Iterate over a the non-zero values from the row and sum to the total of the row
        for (int k = row_start; k < row_end; k++) {
            y[i]+=values[k] *x[col_indices[k]];
        }
    }
}

// =========================================================
// TEST HARNESS
// =========================================================
int main(void) {
    srand(time(NULL));
    
    const int num_iterations = 100;
    int passed_count = 0;

    for (int iter = 0; iter < num_iterations; ++iter) {
        int rows = rand() % 41 + 5;
        int cols = rand() % 41 + 5;
        double density = 0.05 + (rand() / (double) RAND_MAX) * 0.35;
        
        size_t mat_sz = (size_t) rows * cols;

        double* A = calloc(mat_sz, sizeof(double));
        for (size_t i = 0; i < mat_sz; ++i) {
            if (((double) rand() / RAND_MAX) < density) {
                A[i] = ((double) rand() / RAND_MAX) * 20.0 - 10.0;
            }
        }

        double* values = malloc(mat_sz * sizeof(double));
        int* col_indices = malloc(mat_sz * sizeof(int));
        int* row_ptrs = malloc((rows + 1) * sizeof(int));
        double* x = malloc(cols * sizeof(double));
        double* y_user = malloc(rows * sizeof(double));
        double* y_ref = calloc(rows, sizeof(double));
        int out_nnz = 0;

        for (int i = 0; i < cols; ++i) {
            x[i] = ((double) rand() / RAND_MAX) * 20.0 - 10.0;
        }

        for (int i = 0; i < rows; ++i) {
            double sum = 0.0;
            for (int j = 0; j < cols; ++j) {
                sum += A[i * cols + j] * x[j];
            }
            y_ref[i] = sum;
        }

        sparse_multiply(rows, cols, A, x, &out_nnz, values, col_indices, row_ptrs, y_user);

        double max_err = 0.0;
        int passed = 1;
        for (int i = 0; i < rows; ++i) {
            double diff = fabs(y_user[i] - y_ref[i]);
            double tol = 1e-7 + 1e-7 * fabs(y_ref[i]); // Mixed absolute/relative tolerance
            if (diff > tol) {
                max_err = fmax(max_err, diff);
                passed = 0;
            }
        }

        if (passed) {
            passed_count++;
        }

        printf(
            "Iter %2d [%3dx%3d, density=%.2f, nnz=%4d]: %s (Max error: %.2e)\n",
            iter, rows, cols, density, out_nnz, passed ? "PASS" : "FAIL", max_err
        );

        free(A);
        free(values);
        free(col_indices);
        free(row_ptrs);
        free(x);
        free(y_user);
        free(y_ref);
    }

    printf(
        "\n%s (%d/%d iterations passed)\n",
        passed_count == num_iterations ? "All tests passed!" : "Some tests failed.",
        passed_count, num_iterations
    );
           
    return passed_count == num_iterations ? 0 : 1;
}
