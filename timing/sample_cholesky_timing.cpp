#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <chrono>
#include <gtsam/base/Matrix.h>

extern "C" {
#include "cholesky.h"
}

using namespace std;
using namespace Eigen;

float randomFloat() {
    return static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
}

MatrixXf generateSPDMatrix(int n) {
    // MatrixXf A = MatrixXf::Random(n, n);
    // MatrixXf spdMatrix = A.transpose() * A;
    MatrixXf spdMatrix = MatrixXf::Identity(n, n);
    return spdMatrix;
}

MatrixXf generateMNMatrix(int m, int n) {
    // MatrixXf A = MatrixXf::Random(m, n);
    MatrixXf A = MatrixXf::Zero(m, n);
    return A;
}

int main() {
    // seed the random number generator
    srand(static_cast<unsigned int>(time(0)));

    for(int i = 0; i < 1000; i++) {

        int n = rand() % 80 + 24;
        int m = rand() % 200 + n + 1;

        MatrixXf matrix(m, m);

        matrix.block(0, 0, n, n) = generateSPDMatrix(n);
        matrix.block(n, 0, m - n, n) = generateMNMatrix(m - n, n);

        auto start = chrono::high_resolution_clock::now();
        partial_factorization4(matrix.data(), n, m);
        auto end = chrono::high_resolution_clock::now();

        chrono::duration<double, std::nano> ns_double = end - start;

        cout << m << " " << n << " " << ns_double.count() << endl;

    }

}
