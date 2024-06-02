#include <iostream>
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

using namespace std;

int main(int argc, char** argv) {

    int num_timesteps;
    cin >> num_timesteps;

    for(int step = 0; step < num_timesteps; step++) {
        int true_step = step + timestep_start;
        printf("step = %d\n", true_step);

        int num_factors;
        cin >> num_factors;

        // Host data
        vector<int> h_csrRowPtrA = {0};
        vector<int> h_csrColIndA;
        vector<float> h_csrValA;
        vector<int> h_csrRowPtrAT = {0};
        vector<int> h_csrColIndAT;
        vector<float> h_csrValAT;
        vector<float> h_b;
        vector<float> h_x;
        vector<int> h_csrRowPtrD;
        vector<int> h_csrColIndD;
        vector<float> h_csrValD;
        vector<int> h_csrRowPtrH;
        vector<int> h_csrColIndH;
        vector<float> h_csrValH;

        for(int factor = 0; factor < num_factors; factor++) {
            int height, width;
            vector<int> ridx;

            cin >> height >> width;

            for(int i = 0; i < height - 1; i++) {
                int ridx_val;
                cin >> ridx_val;
                ridx.push_back(ridx_val);
            }

            for(int i = 0; i < height - 1; i++) {
                h_csrRowPtrA.push_back(h_csrRowPtrA.back() + width);
                h_csrColIndA.insert(h_csrColIndA.end(), ridx.begin(), ridx.end());
                for(int j = 0; j < width; j++) {
                    if(j == i) {
                        h_csrValA.push_back(1.0);
                    } else {
                        h_csrValA.push_back(0.0);
                    }
                }
            }
        }

        // Transpose A
        h_csrRowPtrAT.resize(remapped_ridx.size() + 1, 0);
        h_csrRowPtrAT[0] = 0;
        for(int i = 0; i < h_csrColIndA.size(); i++) {
            h_csrRowPtrAT[h_csrColIndA[i] + 1]++;
        }

        for(int i = 1; i < h_csrRowPtrAT.size(); i++) {
            h_csrRowPtrAT[i] += h_csrRowPtrAT[i - 1];
        }
        h_csrColIndAT.resize(h_csrColIndA.size());
        h_csrValAT.resize(h_csrValA.size());
        vector<int> h_csrRowPtrAT_tmp = h_csrRowPtrAT;

        for(int i = 0; i < h_csrRowPtrA.size() - 1; i++) {
            for(int j = h_csrRowPtrA[i]; j < h_csrRowPtrA[i + 1]; j++) {
                int col = h_csrColIndA[j];
                int idx = h_csrRowPtrAT_tmp[col]++;
                h_csrColIndAT[idx] = i;
                h_csrValAT[idx] = h_csrValA[j];
            }
        }

	// clock_t start, end;    
	// start = clock();

        // // Convenience variables
        // // A is m x n, AT is n x m, b is m x 1, ATb is n x 1, x is n x 1
        // int m = h_csrRowPtrA.size() - 1;
        // int n = h_csrRowPtrAT.size() - 1;
        // int nnzA = h_csrColIndA.size();
        // float one = 1, zero = 0;

        // // Device memory allocation
        // cudaMalloc(&d_csrRowPtrA, h_csrRowPtrA.size() * sizeof(int));
        // cudaMalloc(&d_csrColIndA, h_csrColIndA.size() * sizeof(int));
        // cudaMalloc(&d_csrValA, h_csrValA.size() * sizeof(float));
        // cudaMalloc(&d_csrRowPtrAT, h_csrRowPtrAT.size() * sizeof(int));
        // cudaMalloc(&d_csrColIndAT, h_csrColIndAT.size() * sizeof(int));
        // cudaMalloc(&d_csrValAT, h_csrValAT.size() * sizeof(float));
        // cudaMalloc(&d_b, h_b.size() * sizeof(float));
        // cudaMalloc(&d_ATb, n * sizeof(float));
        // cudaMalloc(&d_x, n * sizeof(float));

        // // Device memory copy
        // cudaMemcpy(d_csrRowPtrA, h_csrRowPtrA.data(), h_csrRowPtrA.size() * sizeof(int), cudaMemcpyHostToDevice);
        // cudaMemcpy(d_csrColIndA, h_csrColIndA.data(), h_csrColIndA.size() * sizeof(int), cudaMemcpyHostToDevice);
        // cudaMemcpy(d_csrValA, h_csrValA.data(), h_csrValA.size() * sizeof(float), cudaMemcpyHostToDevice);
        // cudaMemcpy(d_csrRowPtrAT, h_csrRowPtrAT.data(), h_csrRowPtrAT.size() * sizeof(int), cudaMemcpyHostToDevice);
        // cudaMemcpy(d_csrColIndAT, h_csrColIndAT.data(), h_csrColIndAT.size() * sizeof(int), cudaMemcpyHostToDevice);
        // cudaMemcpy(d_csrValAT, h_csrValAT.data(), h_csrValAT.size() * sizeof(float), cudaMemcpyHostToDevice);
        // cudaMemcpy(d_b, h_b.data(), h_b.size() * sizeof(float), cudaMemcpyHostToDevice);

        // // printDeviceVals(d_csrRowPtrA, m + 1, "d_csrRowPtrA", "int");
        // // printDeviceVals(d_csrColIndA, nnzA, "d_csrColIndA", "int");
        // // printDeviceVals(d_csrValA, nnzA, "d_csrValA", "float");
        // // printDeviceVals(d_csrRowPtrAT, n + 1, "d_csrRowPtrAT", "int");
        // // printDeviceVals(d_csrColIndAT, nnzA, "d_csrColIndAT", "int");
        // // printDeviceVals(d_csrValAT, nnzA, "d_csrValAT", "float");

        // // Matrix descriptors
        // cusparseSpMatDescr_t descrSpA;
        // cusparseCreateCsr(&descrSpA, m, n, nnzA, d_csrRowPtrA, d_csrColIndA, d_csrValA, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F);
        // cusparseSpMatDescr_t descrSpAT;
        // cusparseCreateCsr(&descrSpAT, n, m, nnzA, d_csrRowPtrAT, d_csrColIndAT, d_csrValAT, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F);
        // cusparseDnVecDescr_t descrDnb;
        // cusparseCreateDnVec(&descrDnb, m, d_b, CUDA_R_32F);
        // cusparseDnVecDescr_t descrDnATb;
        // cusparseCreateDnVec(&descrDnATb, n, d_ATb, CUDA_R_32F);

        // // Compute ATb
        // cusparseSpMV_bufferSize(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
        //                         &one, descrSpAT, descrDnb, &zero, descrDnATb, 
        //                         CUDA_R_32F, CUSPARSE_MV_ALG_DEFAULT, &bufferSize);

        // cudaMalloc(&buffer1, bufferSize);

        // cusparseSpMV(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
        //              &one, descrSpAT, descrDnb, &zero, descrDnATb, CUDA_R_32F, CUSPARSE_CSRMV_ALG1, buffer1);

        // // printDeviceVals(d_b, m, "d_b", "float");
        // // printDeviceVals(d_ATb, n, "ATb", "float");

        // // Compute ATA

        // // Set up zero matrix D
        // int nnzD = 0;
        // h_csrRowPtrD.resize(n + 1, 0);
        // cudaMalloc(&d_csrRowPtrD, h_csrRowPtrD.size() * sizeof(int));
        // cudaMemcpy(d_csrRowPtrD, h_csrRowPtrD.data(), h_csrRowPtrD.size() * sizeof(int), cudaMemcpyHostToDevice);

        // // Matrix descriptors
        // cusparseMatDescr_t descrA, descrAT, descrH, descrD;
        // cusparseCreateMatDescr(&descrA);
        // cusparseCreateMatDescr(&descrAT);
        // cusparseCreateMatDescr(&descrH);
        // cusparseCreateMatDescr(&descrD);

        // // assume matrices A, AT and D are ready
        // int nnzH;
        // csrgemm2Info_t info = NULL;

        // // step 1: create an opaque structure
        // cusparseCreateCsrgemm2Info(&info);

        // cusparseScsrgemm2_bufferSizeExt(cusparseHandle, n, n, m, &one, 
        //                                 descrAT, nnzA, d_csrRowPtrAT, d_csrColIndAT, 
        //                                 descrA, nnzA, d_csrRowPtrA, d_csrColIndA, 
        //                                 &zero, 
        //                                 descrD, nnzD, NULL, NULL,
        //                                 info,
        //                                 &bufferSize);

        // cudaMalloc(&buffer2, bufferSize);

        // // step 3: compute csrRowPtrH
        // cudaMalloc(&d_csrRowPtrH, (n + 1) * sizeof(int));
        // cusparseXcsrgemm2Nnz(cusparseHandle, n, n, m,
        //                      descrAT, nnzA, d_csrRowPtrAT, d_csrColIndAT,
        //                      descrA, nnzA, d_csrRowPtrA, d_csrColIndA,
        //                      descrD, nnzD, d_csrRowPtrD, d_csrColIndD,
        //                      descrH, d_csrRowPtrH, &nnzH,
        //                      info, buffer2);

        // // printDeviceVals(d_csrRowPtrH, n + 1, "d_csrRowPtrH", "int");
        // // printf("nnzH = %d\n", nnzH);

        // // step 4: finish sparsity pattern and value of H
        // cudaMalloc(&d_csrColIndH, nnzH * sizeof(int));
        // cudaMalloc(&d_csrValH, nnzH * sizeof(float));

        // cusparseScsrgemm2(cusparseHandle, n, n, m, &one, 
        //                   descrAT, nnzA, d_csrValAT, d_csrRowPtrAT, d_csrColIndAT, 
        //                   descrA, nnzA, d_csrValA, d_csrRowPtrA, d_csrColIndA, 
        //                   &zero, 
        //                   descrD, nnzD, NULL, d_csrRowPtrD, NULL,
        //                   descrH, d_csrValH, d_csrRowPtrH, d_csrColIndH,
        //                   info, buffer2);

        // // printDeviceVals(d_csrRowPtrH, n + 1, "d_csrRowPtrH", "int");
        // // printDeviceVals(d_csrColIndH, nnzH, "d_csrColIndH", "int");
        // // printDeviceVals(d_csrValH, nnzH, "d_csrValH", "float");

        // // Solve ATA x = ATb
        // int singularity;

        // cusolverSpScsrlsvchol(cusolverSpHandle, n, nnzH, 
        //                       descrH, d_csrValH, d_csrRowPtrH, d_csrColIndH, 
        //                       d_ATb, 0.0, 0, d_x, &singularity);

        // h_x.resize(n);

        // cudaMemcpy(h_x.data(), d_x, n * sizeof(float), cudaMemcpyDeviceToHost);

        // // printf("x: ");
        // // for(int i = 0; i < n; i++) {
        // //     printf("%f ", h_x[i]);
        // // }
        // // printf("\n");

        // cudaFree(d_csrRowPtrA);
        // cudaFree(d_csrColIndA);
        // cudaFree(d_csrValA);
        // cudaFree(d_csrRowPtrAT);
        // cudaFree(d_csrColIndAT);
        // cudaFree(d_csrValAT);
        // cudaFree(d_b);
        // cudaFree(d_ATb);
        // cudaFree(d_x);

        // cudaFree(d_csrRowPtrD);
        // cudaFree(d_csrColIndD);
        // cudaFree(d_csrValD);
        // cudaFree(d_csrRowPtrH);
        // cudaFree(d_csrColIndH);
        // cudaFree(d_csrValH);

        // cudaFree(buffer1);
        // cudaFree(buffer2);

        // cusparseDestroySpMat(descrSpA);
        // cusparseDestroySpMat(descrSpAT);
        // cusparseDestroyDnVec(descrDnb);
        // cusparseDestroyDnVec(descrDnATb);

        // cusparseDestroyMatDescr(descrA);
        // cusparseDestroyMatDescr(descrAT);
        // cusparseDestroyMatDescr(descrH);
        // cusparseDestroyMatDescr(descrD);

        // cusparseDestroyCsrgemm2Info(info);

	// end = clock();
	// double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	// printf("step %d time: %f ms\n", step, cpu_time_used * 1000);
    }

    // cusparseDestroy(cusparseHandle);
    // cusolverSpDestroy(cusolverSpHandle);

    printf("Passed :)\n");
}

