#include <iostream>
#include <Eigen/Core>
#include <unistd.h>
#include <getopt.h>

using namespace std;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> Matrix;

int main(int argc, char** argv) {
    int n = 1000;
    int num_iter = 100;

    // Get experiment setup
    static struct option long_options[] = {
        {"size", required_argument, 0, 'n'},
        {"iter", required_argument, 0, 'i'},
        {0, 0, 0, 0}
    };
    int opt, option_index;
    while((opt = getopt_long(argc, argv, "n:i:", long_options, &option_index)) != -1) {
        switch(opt) {
            case 'n':
                n = atoi(optarg);
                break;
            case 'i':
                num_iter = atoi(optarg);
                break;
        }
    }

    Matrix A(n, n);
    Matrix B(n, n);
    Matrix C(n, n);

    for(int i = 0; i < num_iter; i++) {
        if(i % 2 == 0) {
            C += A * B;
        }
        else {
            C -= A * B;
        }
    }
    
    cout << "Background done" << endl;

    return 0;

}
