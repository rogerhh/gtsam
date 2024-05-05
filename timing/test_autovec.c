#include <stdio.h>

void matmul(int dim_I, int dim_J, int dim_K,
            float* A, float* B, float* C,
            int stride_A, int stride_B, int stride_C) {
  for(size_t i = 0; i < dim_I; i++) {
    for(size_t j = 0; j < dim_J; j++) {
      for(size_t k = 0; k < dim_K; k++) {
        C[i * stride_C + j] = A[i * stride_A + k] * B[k * stride_B + j];
      }
    }
  }
}

int main() {
    int LEN = 128;

    float a[1024] = {1, 2, 3, 4, 5, };
    float b[1024] = {1, 2, 3, 4, 5, };
    float c[1024] = {1, 2, 3, 4, 5, };

    matmul(5, 5, 5, a, b, c, 10, 10, 10);

    for(int i = 0; i < LEN; i++) {
	printf("%f ", a[i]);
    }
    printf("\n");
}
