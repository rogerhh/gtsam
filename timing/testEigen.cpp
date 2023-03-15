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

int main() {
    // Eigen::Matrix<double, C, R> m;
    // Eigen::Matrix<double, C, R> n;
    // // for(int i = 0; i < 3; i++) {
    // //     for(int j = 0; j < 3; j++) {
    // //         m(i, j) = i * j;
    // //     }
    // // }

    // // for(int i = 0; i < 3; i++) {
    // //     for(int j = 0; j < 3; j++) {
    // //         n(i, j) = i - j;
    // //     }
    // // }
    // 
    // long long int sum = 0;
    // m = Eigen::Matrix<double, R, C>::Random();
    // n = Eigen::Matrix<double, C, R>::Random();

    // Eigen::Matrix<double, R, R> k;
    // for(long long int j = 0; j < 10000000; j++) {

    //         m(2, 2) = 0.5;
    //         
    //         // k.noalias() = m.transpose() * n;
    //         k.selfadjointView<Eigen::Upper>().rankUpdate(m.transpose());
    //         // k.selfadjointView<Eigen::Upper>().rankUpdate(m);
    //         // k += m * m.transpose();

    //         // gemm(m.transpose(), m, k);
    // }

    // cout << m << endl << endl << n << endl << endl << k << endl;
    // 
    // // m(1, 2) = 5;

    // // auto m2 = m.selfadjointView<Eigen::Upper>();
    // // // auto m2 = m.triangularView<Eigen::Upper>();
    // // auto m3 = m2 * m;
    // // m = m.selfadjointView<Eigen::Upper>().rankUpdate(d);
    // // auto m4 = m + d;
    // // // m2(2, 1) = 6;
    // // cout << m << endl << endl << d << endl << endl << m4 << endl;
}
