#include <iostream>
#include <gtsam/base/Matrix.h>

#include <gtsam/linear/gemmini_functions.h>

using namespace gtsam;
using namespace std;

using RowMajorMatrix = Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>; 

int main() {
  RowMajorMatrix m(3, 3);
  m << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  RowMajorMatrix n(3, 3);
  n << 9, 8, 7, 6, 5, 4, 3, 2, 1;
  RowMajorMatrix C(4, 4);
  C.setZero();

    for(int i = 0; i < 6; i++) {
        cout << m.data()[i] << endl;
    }


  matmul(3, 3, 3,
         m.data(), n.data(), C.data(),
         3, 3, 4,
         1, 1,
         true, true
	 );

  cout << "m = \n" << m << endl;
  cout << "n = \n" << n << endl;
  cout << "C = \n" << C << endl;


}
