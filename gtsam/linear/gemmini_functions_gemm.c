/**
* @file    gemmini_functions_gemm.cpp
* @brief   All the different linear algebra routines used in CholeskyEliminationTree. Isolated out to support gemmini
* @author  Roger Hsiao (rogerhh)
* @date    Aug. 11, 2023
*/

#include <gtsam/linear/gemmini_functions.h>
#include <assert.h>

// Perform C += A_scale_factor * A^(transpose_A * T) * B_scale_factor * B^(transpose_B * T)
// Assume A, B, C are stored in row-major order
void matmul(
  size_t dim_I, size_t dim_J, size_t dim_K, 
  const elem_t* A, const elem_t* B, elem_t* C,
  size_t stride_A, size_t stride_B, size_t stride_C,
  scale_t A_scale_factor, scale_t B_scale_factor,
  bool transpose_A, bool transpose_B) {

  tiled_matmul_auto(dim_I, dim_J, dim_K,
                    A, B, C, C,
                    stride_A, stride_B, stride_C, stride_C,
                    A_scale_factor, B_scale_factor, 1,
                    0, 1, 0,
                    false,
                    transpose_A, transpose_B,
                    false, false,
                    0,
                    CPU);

}

// Perform C = A_scale_factor * A^(transpose_A * T) * B_scale_factor * B^(transpose_B * T) + D_scale_factor * D
// Assume A, B, C, D are stored in row-major order
void matmul2(
  size_t dim_I, size_t dim_J, size_t dim_K, 
  const elem_t* A, const elem_t* B, 
  const elem_t* D, elem_t* C,
  size_t stride_A, size_t stride_B, size_t stride_D, size_t stride_C,
  scale_t A_scale_factor, scale_t B_scale_factor, scale_t D_scale_factor,
  bool transpose_A, bool transpose_B) {

  tiled_matmul_auto(dim_I, dim_J, dim_K,
                    A, B, D, C,
                    stride_A, stride_B, stride_D, stride_C,
                    A_scale_factor, B_scale_factor, D_scale_factor,
                    0, 1, 0,
                    false,
                    transpose_A, transpose_B,
                    false, false,
                    0,
                    CPU);

}

// Perform C += A_scale_factor * A^(transpose_A * T) * A^(transpose_A * T)
// Assume A, C are stored in row-major order
// Inputting B as a parameter because I don't want to change the function call lol
// B is required to be the same as A
void syrk(
  size_t dim_I, size_t dim_J, size_t dim_K, 
  const elem_t* A, const elem_t* B, elem_t* C,
  size_t stride_A, size_t stride_B, size_t stride_C,
  scale_t A_scale_factor, scale_t B_scale_factor,
  bool transpose_A, bool transpose_B) {
  size_t stride_Ai, stride_Ak, stride_Bk, stride_Bj, stride_Ci = stride_C, stride_Cj = 1;

  tiled_matmul_auto(dim_I, dim_J, dim_K,
                    A, B, C, C,
                    stride_A, stride_B, stride_C, stride_C,
                    A_scale_factor, B_scale_factor, 1,
                    0, 1, 0,
                    false,
                    transpose_A, transpose_B,
                    false, false,
                    0,
                    CPU);

  /*
  if(transpose_A) {
    stride_Ai = 1;
    stride_Ak = stride_A;
  }
  else {
    stride_Ai = stride_A;
    stride_Ak = 1;
  }
  if(transpose_B) {
    stride_Bk = 1;
    stride_Bj = stride_B;
  }
  else {
    stride_Bk = stride_B;
    stride_Bj = 1;
  }
  const elem_t* Ai = A;
  elem_t* Ci = C;
  for(size_t i = 0; i < dim_I; i++) {
    const elem_t* Bj = B + i * stride_Bj;
    elem_t* Cij = Ci + i * stride_Cj;
    for(size_t j = i; j < dim_J; j++) {
      const elem_t* Aik = Ai;
      const elem_t* Bkj = Bj;
      for(size_t k = 0; k < dim_K; k++) {
        *Cij += A_scale_factor * (*Aik) * B_scale_factor * (*Bkj);
        Aik += stride_Ak;
        Bkj += stride_Bk;
      }
      Bj += stride_Bj;
      Cij += stride_Cj;
    }
    Ai += stride_Ai;
    Ci += stride_Ci;
  }
  */
}

// Perform y += A_scale_factor A * x
void gemv(
  size_t dim_I, size_t dim_K,
  const elem_t* A, const elem_t* x, elem_t* y,
  size_t stride_A,
  scale_t A_scale_factor) {
  size_t stride_Ai = stride_A, stride_Ak = 1, stride_xk = 1, stride_yi = 1;

  const elem_t* Ai = A;
  elem_t* yi = y;
  for(size_t i = 0; i < dim_I; i++) {
    const elem_t* Aik = Ai;
    const elem_t* xk = x;
    for(size_t k = 0; k < dim_K; k++) {
      *yi += A_scale_factor * (*Aik) * (*xk);
      Aik += stride_Ak;
      xk += stride_xk;
    }
    Ai += stride_Ai;
    yi += stride_yi;
  }
}

// Perform y = A_scale_factor A * x_scale_factor x + z_scale_factor z
// A is of dim_I x dim_K
// x is a vector of dim_K x 1, y, z are vectors of dim_I x 1
void gemv2(
  size_t dim_I, size_t dim_K,
  const elem_t* A, const elem_t* x, 
  const elem_t* z, elem_t* y,
  size_t stride_A,
  scale_t A_scale_factor, scale_t x_scale_factor, scale_t z_scale_factor,
  bool transpose_A, bool transpose_x) {

  tiled_matmul_auto(dim_I, 1, dim_K,
                    A, x, z, y,
                    stride_A, 1, 1, 1,
                    A_scale_factor, x_scale_factor, z_scale_factor,
                    0, 1, 0,
                    false,
                    transpose_A, transpose_x,
                    false, false,
                    0,
                    CPU);
}

// // Do C += scale * A * B. A, B, C must be preallocated
// void matmul(
//     size_t A_rows, size_t A_cols, 
//     size_t B_rows, size_t B_cols,
//     GEMMINI_TYPE scale,
//     const GEMMINI_TYPE* A, const GEMMINI_TYPE* B, GEMMINI_TYPE* C) {
//   assert(A_cols == B_rows);
//   for(size_t i = 0; i < A_rows; i++) {
//     for(size_t j = 0; j < B_rows; j++) {
//       for(size_t k = 0; k < A_cols; k++) {
//         C[j * A_rows + i] += scale * A[k * A_rows + i] * B[j * B_rows + k];
//       }
//     }
//   }
// }

// Do C += scale * A.T * B. A, B, C must be preallocated
void transposed_matmul(
    size_t A_rows, size_t A_cols, 
    size_t B_rows, size_t B_cols,
    GEMMINI_TYPE scale,
    const GEMMINI_TYPE* A, const GEMMINI_TYPE* B, GEMMINI_TYPE* C) {
  assert(A_rows == B_rows);
  for(size_t i = 0; i < A_cols; i++) {
    for(size_t j = 0; j < B_rows; j++) {
      for(size_t k = 0; k < A_rows; k++) {
        C[j * A_cols + i] += scale * A[k * A_cols + i] * B[j * B_rows + k];
      }
    }
  }
}

// // Do C += scale * A.T * A. A, C must be preallocated
// void syrk(
//   size_t A_rows, size_t A_cols,
//   GEMMINI_TYPE scale,
//   const GEMMINI_TYPE* A, GEMMINI_TYPE* C) {
//   for(size_t j = 0; j < A_cols; j++) {
//     GEMMINI_TYPE* C_start = C + j * A_cols;
//     const GEMMINI_TYPE* A_start = A + j * A_rows;
//     for(size_t i = j; i < A_cols; i++) {
//       const GEMMINI_TYPE* AT_start = A + i * A_rows;
//       for(size_t k = 0; k < A_rows; k++) {
//         *(C_start + i) += scale * (*(AT_start + k)) * (*(A_start + k));
//       }
//     }
//   }
// }

// // Do y += scale * Ax
// void gemv(
//   size_t A_rows, size_t A_cols,
//   GEMMINI_TYPE scale,
//   const GEMMINI_TYPE* A, 
//   const GEMMINI_TYPE* x,
//   GEMMINI_TYPE* y) {
// 
//   for(size_t j = 0; j < A_cols; j++) {
//     const GEMMINI_TYPE* A_start = A + j * A_rows;
//     for(size_t i = 0; i < A_rows; i++) {
//       y[i] += scale * A_start[i] * x[j];
//     }
//   }
// 
// }

// Do y += scale * A.Tx
void transpose_gemv(
  size_t A_rows, size_t A_cols,
  GEMMINI_TYPE scale,
  const GEMMINI_TYPE* A, 
  const GEMMINI_TYPE* x, 
  GEMMINI_TYPE* y) {

  for(size_t j = 0; j < A_cols; j++) {
    const GEMMINI_TYPE* A_start = A + j * A_rows;
    for(size_t i = 0; i < A_rows; i++) {
      y[j] += scale * A_start[i] * x[i];
    }
  }
}

// Add the block of A(Ar, Ac, h, w) into the block of M(Mr, Mc, h, w)
void add_block(
  const GEMMINI_TYPE* A, GEMMINI_TYPE* M,
  size_t stride_A, size_t stride_M,
  size_t Ar, size_t Ac,
  size_t Mr, size_t Mc,
  size_t h, size_t w
  ) {
  const GEMMINI_TYPE* A_col_start = A + Ac * stride_A + Ar;
  GEMMINI_TYPE* M_col_start = M + Mc * stride_M + Mr;
  for(size_t j = 0; j < w; j++) {
    const GEMMINI_TYPE* A_ptr = A_col_start;
    GEMMINI_TYPE* M_ptr = M_col_start;
    for(size_t i = 0; i < h; i++) {
      *M_ptr += *A_ptr;
      A_ptr++;
      M_ptr++;
    }
    A_col_start += stride_A;
    M_col_start += stride_M;
  }
}

// Add the block of A(Ar, Ac, h, w).T into the block of M(Mr, Mc, w, h)
void add_block_transpose(
  const GEMMINI_TYPE* A, GEMMINI_TYPE* M,
  size_t stride_A, size_t stride_M,
  size_t Ar, size_t Ac,
  size_t Mr, size_t Mc,
  size_t h, size_t w
  ) {
  const GEMMINI_TYPE* A_col_start = A + Ac * stride_A + Ar;
  GEMMINI_TYPE* M_col_start = M + Mc * stride_M + Mr;
  for(size_t j = 0; j < h; j++) {
    const GEMMINI_TYPE* A_ptr = A_col_start;
    GEMMINI_TYPE* M_ptr = M_col_start;
    for(size_t i = 0; i < w; i++) {
      *M_ptr += *A_ptr;
      A_ptr += stride_A;
      M_ptr++;
    }
    A_col_start++;
    M_col_start += stride_M;
  }
}
