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

void printDeviceVals(float* d_vals, int n, const string& name, const string& type) {
    if(type == "int") {
        int* h_vals = (int*)malloc(n * sizeof(int));
        cudaMemcpy(h_vals, d_vals, n * sizeof(int), cudaMemcpyDeviceToHost);
        printf("%s: ", name.c_str());
        for(int i = 0; i < n; i++) {
            printf("%d ", h_vals[i]);
        }
        printf("\n");
        free(h_vals);
    }
    else if(type == "float") {
        float* h_vals = (float*)malloc(n * sizeof(float));
        cudaMemcpy(h_vals, d_vals, n * sizeof(float), cudaMemcpyDeviceToHost);
        printf("%s: ", name.c_str());
        for(int i = 0; i < n; i++) {
            printf("%f ", h_vals[i]);
        }
        printf("\n");
        free(h_vals);
    }
}

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

        // Host data
        vector<int> h_csrRowPtrA = {0};
        vector<int> h_csrColIndA;
        vector<float> h_csrValA;
        vector<int> h_csrRowPtrAT = {0};
        vector<int> h_csrColIndAT;
        vector<float> h_csrValAT;
        vector<float> h_b;

        // Device data
        int* d_csrRowPtrA;
        int* d_csrColIndA;
        float* d_csrValA;
        int* d_csrRowPtrAT;
        int* d_csrColIndAT;
        float* d_csrValAT;
        float* d_b;
        float* d_ATb;
        float* d_x;
        size_t bufferSize;
        void* buffer = NULL;

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
                    ridx_set.insert(node_ridx[node][ridx[ih]]);
                }

                for(int j = 0; j < width; j++) {
                    h_b.push_back(1.0f);
                    h_csrRowPtrA.push_back(h_csrRowPtrA.back() + height - 1);
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

        // Transpose A
        h_csrRowPtrAT.resize(remapped_ridx.size() + 1, 0);
        for(int i = 0; i < h_csrColIndA.size(); i++) {
            h_csrRowPtrAT[h_csrColIndA[i] + 1]++;
        }
        for(int i = 1; i < h_csrRowPtrAT.size(); i++) {
            h_csrRowPtrAT[i] += h_csrRowPtrAT[i - 1];
        }
        h_csrColIndAT.resize(h_csrColIndA.size());
        h_csrValAT.resize(h_csrValA.size());
        for(int i = 0; i < h_csrRowPtrA.size() - 1; i++) {
            for(int j = h_csrRowPtrA[i]; j < h_csrRowPtrA[i + 1]; j++) {
                int col = h_csrColIndA[j];
                int idx = h_csrRowPtrAT[col]++;
                h_csrColIndAT[idx] = i;
                h_csrValAT[idx] = h_csrValA[j];
            }
        }

        // Convenience variables
        // A is m x n, AT is n x m, b is m x 1, ATb is n x 1, x is n x 1
        int m = h_csrRowPtrA.size() - 1;
        int n = h_csrRowPtrAT.size() - 1;
        int nnzA = h_csrColIndA.size();
        int one = 1, zero = 0;

        // Device memory allocation
        cudaMalloc(&d_csrRowPtrA, h_csrRowPtrA.size() * sizeof(int));
        cudaMalloc(&d_csrColIndA, h_csrColIndA.size() * sizeof(int));
        cudaMalloc(&d_csrValA, h_csrValA.size() * sizeof(float));
        cudaMalloc(&d_csrRowPtrAT, h_csrRowPtrAT.size() * sizeof(int));
        cudaMalloc(&d_csrColIndAT, h_csrColIndAT.size() * sizeof(int));
        cudaMalloc(&d_csrValAT, h_csrValAT.size() * sizeof(float));
        cudaMalloc(&d_b, h_b.size() * sizeof(float));
        cudaMalloc(&d_ATb, n * sizeof(float));
        cudaMalloc(&d_x, n * sizeof(float));

        // Device memory copy
        cudaMemcpy(d_csrRowPtrA, h_csrRowPtrA.data(), h_csrRowPtrA.size() * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_csrColIndA, h_csrColIndA.data(), h_csrColIndA.size() * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_csrValA, h_csrValA.data(), h_csrValA.size() * sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_csrRowPtrAT, h_csrRowPtrAT.data(), h_csrRowPtrAT.size() * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_csrColIndAT, h_csrColIndAT.data(), h_csrColIndAT.size() * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_csrValAT, h_csrValAT.data(), h_csrValAT.size() * sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_b, h_b.data(), h_b.size() * sizeof(float), cudaMemcpyHostToDevice);

        // Matrix descriptors
        cusparseSpMatDescr_t descrA;
        cusparseCreateCsr(&descrA, m, n, nnzA, d_csrRowPtrA, d_csrColIndA, d_csrValA, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F);
        cusparseSpMatDescr_t descrAT;
        cusparseCreateCsr(&descrAT, n, m, nnzA, d_csrRowPtrAT, d_csrColIndAT, d_csrValAT, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F);
        cusparseDnVecDescr_t descrb;
        cusparseCreateDnVec(&descrb, m, d_b, CUDA_R_32F);
        cusparseDnVecDescr_t descrATb;
        cusparseCreateDnVec(&descrATb, n, d_ATb, CUDA_R_32F);

        // Compute ATb
        cusparseSpMV_bufferSize(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                &one, descrAT, descrb, &zero, descrATb, 
                                CUDA_R_32F, CUSPARSE_CSRMV_ALG1, &bufferSize);

        printf("bufferSize = %lu\n", bufferSize);

        cudaMalloc(&buffer, bufferSize);

        cusparseSpMV(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                     &one, descrAT, descrb, &zero, descrATb, CUDA_R_32F, CUSPARSE_CSRMV_ALG1, buffer);

        printDeviceVals(d_ATb, n, "ATb", "float");

	end = clock();
	double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("step %d time: %f ms\n", step, cpu_time_used * 1000);

        exit(1);
    }

    printf("Passed :)\n");
}

