#include <stdio.h>
#include <gtsam/base/Matrix.h>

int main() {
    int LEN = 128;

    double a[1024] = {0};


    for(int i = 0; i < LEN; i++) {
	a[i] = i;
    }

    for(int i = 0; i < LEN; i++) {
	a[i] *= 2;
    }

    for(int i = 0; i < LEN; i++) {
	printf("%f ", a[i]);
    }
    printf("\n");
}
