#include <iostream>
#include <gtsam/base/Matrix.h>

extern "C" {
#include <gtsam/linear/gemmini_functions.h>
}

using namespace gtsam;
using namespace std;

template<typename MATRIX1, typename MATRIX2, typename MATRIX3>
void call_matmul(
  size_t dim_I, size_t dim_J, size_t dim_K, 
  MATRIX1& A, MATRIX2& B, MATRIX3& C,
  size_t stride_A, size_t stride_B, size_t stride_C,
  scale_t A_scale_factor, scale_t B_scale_factor,
  bool transpose_A, bool transpose_B) {

  matmul(dim_I, dim_J, dim_K, 
         &A(0, 0), &B(0, 0), &C(0, 0), 
         stride_A, stride_B, stride_C, 
         A_scale_factor, B_scale_factor, 
         transpose_A, transpose_B);
}

using RowMajorMatrix = Eigen::Matrix<GEMMINI_TYPE, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>; 

int main() {
  RowMajorMatrix m(3, 3);
  m << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  RowMajorMatrix n(3, 3);
  n << 2, 4, 6, 3, 5, 7, 1, 6, 2;
  RowMajorMatrix C(3, 3);
  C.setZero();

  matmul(2, 2, 2, &m(0, 0), &n(0, 0), &C(0, 0), 3, 3, 3, 1, 1, false, false);

  cout << "m = \n" << m << endl << endl;
  cout << "n = \n" << n << endl << endl;
  cout << "C = \n" << C << endl;

  matmul(3, 3, 2, &m(0, 0), &n(0, 0), &C(0, 0), 3, 3, 3, 1, 1, false, false);

  cout << "m = \n" << m << endl << endl;
  cout << "n = \n" << n << endl << endl;
  cout << "C = \n" << C << endl;

  auto m_map = Eigen::Map<RowMajorMatrix>(&m(0, 0), 3, 3);

  call_matmul(2, 2, 2, m_map, n, C, 3, 3, 3, 1, 1, false, false);
  cout << "m = \n" << m << endl << endl;
  cout << "n = \n" << n << endl << endl;
  cout << "C = \n" << C << endl;

  cout << "Z_3x3: \n" << Z_3x3 << endl;

}
