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
#include <set>
#include <map>

#include "baremetal_tests/incremental_sphere2500_steps-2-200_period-25/incremental_dataset.h"

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

        vector<float> h_b;

        set<int> ridx_set;
        map<int, int> remapped_ridx;

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

                for(int ih = 0; ih < height - 1; ih++) {
                    ridx_set.insert(ridx[ih]);
                }

                for(int j = 0; j < width; j++) {
                    h_b.push_back(0.0f);
                    h_csrRowPtrA.push_back(h_csrRowPtrA.back() + height);
                    for(int ih = 0; ih < height - 1; ih++) {
                        printf("%d %d\n", ih, ridx[ih]);
                        h_csrColIndA.push_back(node_ridx[node][ridx[ih]]);
                        h_csrValA.push_back(1.0f);
                    }
                }
            }
        }

        int count = 0;
        for(int ridx : ridx_set) {
            remapped_ridx[ridx] = count++;
        }

        for(int i = 0; i < h_csrColIndA.size(); i++) {
            h_csrColIndA[i] = remapped_ridx[h_csrColIndA[i]];
        }

        printf("h_csrRowPtrA: ");
        for(int i = 0; i < h_csrRowPtrA.size(); i++) {
            printf("%d ", h_csrRowPtrA[i]);
        }
        printf("\n");
        printf("h_csrColIndA: ");
        for(int i = 0; i < h_csrColIndA.size(); i++) {
            printf("%d ", h_csrColIndA[i]);
        }
        printf("\n");
        printf("h_csrValA: ");
        for(int i = 0; i < h_csrValA.size(); i++) {
            printf("%f ", h_csrValA[i]);
        }
        printf("\n");

	end = clock();
	double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("step %d time: %f ms\n", step, cpu_time_used * 1000);

        exit(1);
    }

    printf("Passed :)\n");
}

