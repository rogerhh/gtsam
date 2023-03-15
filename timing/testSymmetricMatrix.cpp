#include <gtsam/base/SparseSymmetricBlockMatrix.h>

using namespace std;
using namespace gtsam;

int main() {
    SparseSymmetricBlockMatrix hessian;
    hessian.addColumn(0, 6);
    hessian.addColumn(1, 3);
    hessian.preallocateBlock(0, 1);

    hessian.resolveAllocate();

    hessian.print(cout);

    cout << hessian.colBlockRange(0, 0, hessian.colHeight(0)) << endl;

}
