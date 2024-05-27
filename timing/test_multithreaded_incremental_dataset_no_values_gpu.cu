#include <cuda_runtime.h>
#include <cusparse_v2.h>
#include <cusolverSp.h>

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdint.h>
#include <time.h>
#include <cstring>
#include <vector>

#include "baremetal_tests/incremental_sphere2500_steps-2-2000_period-25/incremental_dataset.h"

#include "memory.h"
#include "cholesky.h"

using namespace std;

int main(int argc, char** argv) {
    // Initialize cuSPARSE and cuSOLVER
    cusparseHandle_t cusparseHandle = NULL;
    cusparseCreate(&cusparseHandle);
    cusparseSetPointerMode(cusparseHandle, CUSPARSE_POINTER_MODE_HOST);

    cusolverSpHandle_t cusolverSpHandle = NULL;
    cusolverSpCreate(&cusolverSpHandle);

    for(int step = 0; step < num_timesteps; step++) {
	clock_t start, end;    
	start = clock();
        int true_step = step + timestep_start;
        printf("step = %d\n", true_step);

        int nnodes = step_nnodes[step];

        bool* node_marked = step_node_marked[step];
        bool* node_fixed = step_node_fixed[step];

        int** node_ridx = step_node_ridx[step];

        int* node_num_factors = step_node_num_factors[step];
        int** node_factor_height = step_node_factor_height[step];
        int** node_factor_width = step_node_factor_width[step];
        int*** node_factor_ridx = step_node_factor_ridx[step];

        // Construct A matrix from factors
        vector<int> h_csrRowPtrA = {0};
        vector<int> h_csrColIndA;
        vector<float> h_csrValA;

        for(int node = 0; node < nnodes - 1; node++) {
            bool marked = node_marked[node];
            bool fixed = node_fixed[node];

            if(!marked && !fixed) { continue; }

            int num_factors = node_num_factors[node];
            int* factor_height = node_factor_height[node];
            int* factor_width = node_factor_width[node];
            int** factor_ridx = node_factor_ridx[node];

            for(int i = 0; i < num_factors; i++) {
                int height = factor_height[i];
                int width = factor_width[i];
                int* ridx = factor_ridx[i];

                for(int ih = 0; ih < height; ih++) {
                    h_csrRowPtrA.push_back(h_csrRowPtrA.back() + width);
                    for(int j = 0; j < width; j++) {
                        h_csrColIndA.push_back(node_ridx[node][ridx[j]]);
                        h_csrValA.push_back(1.0f);
                    }
                }
            }
        }


	end = clock();
	double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("step %d time: %f ms\n", step, cpu_time_used * 1000);
    }

    printf("Passed :)\n");
}

