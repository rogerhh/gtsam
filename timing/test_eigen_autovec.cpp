#include <iostream>
#include <gtsam/base/Matrix.h>
#include <gtsam/linear/gemmini_functions.h>

using namespace gtsam;
using namespace std;

using ColMajorMatrix = Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>; 

int main() {
  ColMajorMatrix m(3, 4);
  m << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;
  ColMajorMatrix n(4, 3);
  n << 2, 4, 6, 3, 5, 7, 1, 6, 2, 2, 4, 1;
  ColMajorMatrix C(3, 3);
  C.setZero();

  matmul(3, 3, 4,
         
	 )
  const elem_t* A, const elem_t* B, elem_t* C,
  size_t stride_A, size_t stride_B, size_t stride_C,
  scale_t A_scale_factor, scale_t B_scale_factor,
  bool transpose_A, bool transpose_B) {



}
