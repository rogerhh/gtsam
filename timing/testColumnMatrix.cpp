#include <gtsam/base/SparseColumnBlockMatrix.h>

using namespace std;
using namespace gtsam;

int main() {
    SparseColumnBlockMatrix cols(1, 6);
    cols.resolveAllocate();
    cols.preallocateBlock(2, 3, true);
    cols.preallocateBlock(3, 8, true);
    cols.resolveAllocate();
    cols.block(3).transpose()(1, 3) = 5.5;
    cols.print(cout);
    cols.setZero();
    cols.print(cout);
    cols.checkInvariant();
}
