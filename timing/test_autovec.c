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
  float* Acol = A;
  float* Ccol = C;
  float AB_scale_factor = 1.0;
  for(size_t i = 0; i < dim_I; i++) {
    for(size_t k = 0; k < dim_K; k++) {
      float* Bcol = B + k;
      for(size_t j = 0; j < dim_J; j++) {
        Ccol[j] += AB_scale_factor * Acol[k] * (*Bcol);
        Bcol += stride_B;
      }
    }
    Acol += stride_A;
    Ccol += stride_C;
  }
}

void matmul3(int dim_I, int dim_J, int dim_K,
            float* A, float* B, float* C,
            int stride_A, int stride_B, int stride_C) {
  float* Acol = A;
  float* Bcol = B;
  float AB_scale_factor = 1.0;
  for(size_t k = 0; k < dim_K; k++) {
    float* Ccol = C;
    for(size_t i = 0; i < dim_I; i++) {
      for(size_t j = 0; j < dim_J; j++) {
        Ccol[j] += AB_scale_factor * Acol[i] * Bcol[j];
      }
      Ccol += stride_C;
    }
    Acol += stride_A;
    Bcol += stride_B;
  }
}

void matmul4(int dim_I, int dim_J, int dim_K,
            float* A, float* B, float* C,
            int stride_A, int stride_B, int stride_C) {
  float* Acol = A;
  float AB_scale_factor = 1.0;
  for(size_t k = 0; k < dim_K; k++) {
    float* Bcol = B + k;
    float* Ccol = C;
    for(size_t i = 0; i < dim_I; i++) {
      for(size_t j = 0; j < dim_J; j++) {
        Ccol[j] += AB_scale_factor * Acol[i] * (*Bcol);
        Bcol += stride_B;
      }
      Ccol += stride_C;
    }
    Acol += stride_A;
  }
} 

void gemv1(int dim_I, int dim_K,
            float* A, float* x, 
            float* z, float* y,
            float z_scale_factor, float Ax_scale_factor,
            int stride_A) {
  for(size_t k = 0; k < dim_K; k++) {
    float* Acol = A + k;
    for(size_t i = 0; i < dim_I; i++) {
      y[i] += Ax_scale_factor * (*Acol) * x[k] + z_scale_factor * z[i];
      Acol += stride_A;
    }
  }
}

void gemv2(int dim_I, int dim_K,
            float* A, float* x, 
            float* z, float* y,
            float z_scale_factor, float Ax_scale_factor,
            int stride_A) {

  float* Acol = A;
  for(size_t k = 0; k < dim_K; k++) {
    for(size_t i = 0; i < dim_I; i++) {
      y[i] += Ax_scale_factor * Acol[i] * x[k] + z_scale_factor * z[i];
    }
    Acol += stride_A;
  }
}

int main() {
    int LEN = 128;

    float a[1024] = {1, 2, 3, 4, 5, };
    float b[1024] = {1, 2, 3, 4, 5, };
    float c[1024] = {1, 2, 3, 4, 5, };

    matmul1(5, 5, 5, a, b, c, 10, 10, 10);
    matmul2(5, 5, 5, a, b, c, 10, 10, 10);
    matmul3(5, 5, 5, a, b, c, 10, 10, 10);
    matmul4(5, 5, 5, a, b, c, 10, 10, 10);

    gemv1(5, 3, a, b, c, c, 1, 1, 10);
    gemv2(5, 3, a, b, c, c, 1, 1, 10);

    for(int i = 0; i < LEN; i++) {
	printf("%f ", a[i]);
    }
    printf("\n");
}
