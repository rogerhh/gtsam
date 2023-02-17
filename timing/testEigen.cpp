#define EIGEN_USE_MKL_ALL
#include <Eigen/Core>  
#include <iostream>
#include <chrono>


typedef Eigen::MatrixXd Matrix;
using namespace std;

const int R = 6, C = 6;

inline void gemm(const Eigen::Matrix<double, R, C>& m, 
                 const Eigen::Matrix<double, C, R>& n, 
                 Eigen::Matrix<double, R, R>& k) {
    for(int i = 0; i < R; i++) {
        for(int j = 0; j < C; j++) {
            // k(i, j) = 0;
            for(int q = 0; q < C; q++) {
                k(i, j) = k(i, j) + m(i, q) * n(q, j);
            }
        }
    }
}

template<typename Xpr>
void func(Eigen::MatrixBase<Xpr>& b) {
    cout << b(2, 1) << endl << endl;
}

int main() {
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m;
    Eigen::Matrix<double, C, R> n;
    // m.conservativeResize(R, C);
    m = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>(R, C);
    for(int i = 1; i <= R; i++) {
        for(int j = 1; j <= C; j++) {
            m(i-1, j-1) = i * 10 + j;
        }
    }
    // Eigen::Block<Matrix> b = m.block(0, 0, 4, 3);
    Eigen::Block<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, -1, -1, false> b = m.block(0, 0, 4, 3);
    // auto b1 = m.block(0, 0, 4, 3);
    // b.test();

    // b(1, 2) = 0;

    // func(b);
    // func(b1);

    cout << m << endl << endl << b << endl;

}
