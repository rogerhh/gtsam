#include <iostream>
#include <gtsam/base/Matrix.h>

#include <gtsam/linear/gemmini_functions.h>

using namespace gtsam;
using namespace std;

using ColMajorMatrix = Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>; 

int main() {
  ColMajorMatrix m(3, 2);
  m << 1, 2, 3, 4, 5, 6;
  ColMajorMatrix n(3, 2);
  n << 6, 5, 4, 3, 2, 1;
  ColMajorMatrix C(4, 4);
  C.setZero();


  matmul(2, 2, 3,
         m.data(), n.data(), C.data(),
         3, 3, 4,
         1, 1,
         false, true
	 );

  cout << "m = \n" << m << endl;
  cout << "n = \n" << n << endl;
  cout << "C = \n" << C << endl;


}
