#include <stdio.h>

void matmul1(int dim_I, int dim_J, int dim_K,
            float* A, float* B, float* C,
            int stride_A, int stride_B, int stride_C) {
  float* Acol = A;
  float* Ccol = C;
  float AB_scale_factor = 1.0;
  for(size_t i = 0; i < dim_I; i++) {
    float* Bcol = B;
    for(size_t k = 0; k < dim_K; k++) {
      for(size_t j = 0; j < dim_J; j++) {
        Ccol[j] += AB_scale_factor * Acol[k] * B[j];
      }
      Bcol += stride_B;
    }
    Acol += stride_A;
    Ccol += stride_C;
  }
}

void matmul2(int dim_I, int dim_J, int dim_K,
            float* A, float* B, float* C,
            int stride_A, int stride_B, int stride_C) {
  float AB_scale_factor = 1.0;
  for(size_t i = 0; i < dim_I; i++) {
    for(size_t j = 0; j < dim_J; j++) {
      for(size_t k = 0; k < dim_K; k++) {
          C[i * stride_C + j] += AB_scale_factor * A[i * stride_A + k] * B[j * stride_B + k];
      }
    }
  }
}

int main() {
    int LEN = 128;

    float a[1024] = {1, 2, 3, 4, 5, };
    float b[1024] = {1, 2, 3, 4, 5, };
    float c[1024] = {1, 2, 3, 4, 5, };

    matmul1(5, 5, 5, a, b, c, 10, 10, 10);
    matmul2(5, 5, 5, a, b, c, 10, 10, 10);

    for(int i = 0; i < LEN; i++) {
	printf("%f ", a[i]);
    }
    printf("\n");
}
