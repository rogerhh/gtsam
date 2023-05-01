/**
* @file    Workspace.h
* @brief   Work space memory for multifrontal elimination
* @author  Roger Hsiao (rogerhh)
* @date    Apr. 12, 2023
*/

#include <Eigen/Dense>
#include <Eigen/Core>
#include <cassert>
#include <vector>
#include <iostream>
#include <utility>
#include <gtsam/base/LowerTriangularColumnMatrix.h>
#include <tuple>

namespace gtsam {

class Workspace {
public:
    using BlockIndexVector = LowerTriangularColumnMatrix::BlockIndexVector;
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> ColMajorMatrix;

    struct BlockInfo {

        BlockInfo(double* ptr_in, size_t r_in, size_t c_in, 
                    size_t cliqueSize_in,
                    BlockIndexVector blockIndices_in)
        : ptr(ptr_in), r(r_in), c(c_in), cliqueSize(cliqueSize_in) {
            blockIndices = std::move(blockIndices_in);
        }

        double* ptr = nullptr;
        size_t r = 0;
        size_t c = 0;
        
        size_t cliqueSize = 0;
        BlockIndexVector blockIndices;
    };

private:

    double* data = nullptr;
    double* stack_ptr = nullptr;
    size_t size = 0;
    std::vector<BlockInfo> blockInfos;

    std::vector<double*> stack_ptrs;
    std::vector<std::vector<Eigen::Map<ColMajorMatrix>>> stack;

public:


    Workspace() {}

    ~Workspace() {
        deallocate();
    }

    void deallocate() {
        free(data);
        data = nullptr;
        size = 0;
    }

    // allocate memory for underlying memory
    // This should only be run once after symbolic elim
    void allocate(size_t size_in) {
        assert(blockInfos.empty());
        data = (double*) realloc(data, size_in * sizeof(double));
        stack_ptr = data;
        size = size_in;
    }

    size_t push_matrices(const std::vector<std::pair<size_t, size_t>>& matrixSizes) {
        stack.push_back(std::vector<Eigen::Map<ColMajorMatrix>>());
        stack_ptrs.push_back(stack_ptr);
        for(const auto& p : matrixSizes) {
            size_t r, c;
            std::tie(r, c) = p;
            assert(stack_ptr + r * c <= data + size);
            stack.back().push_back(Eigen::Map<ColMajorMatrix>(stack_ptr, r, c));
            stack.back().back().setZero();
            stack_ptr += r * c;
        }
        return stack.size() - 1;
    }

    const std::vector<Eigen::Map<ColMajorMatrix>>& get_matrices(size_t index) {
        return stack.at(index);
    }

    void pop_matrices() {
        stack_ptr = stack_ptrs.back();
        stack.pop_back();
        stack_ptrs.pop_back();
    }

    // Returns index of matrix on the stack
    size_t push_matrix(size_t r, size_t c, size_t cliqueSize, BlockIndexVector& blockIndices) {
        blockInfos.push_back(BlockInfo(stack_ptr, r, c, cliqueSize, blockIndices));
        stack_ptr += r * c;
        assert(stack_ptr <= data + size);
        return blockInfos.size() - 1;
    }

    // pop the last matrix on the stack
    void pop_matrix() {
        stack_ptr = blockInfos.back().ptr;
        blockInfos.pop_back();
    }

    // Returns matrix at the index
    Eigen::Map<ColMajorMatrix> matrix(size_t index) {
        BlockInfo b = blockInfos[index];
        return Eigen::Map<ColMajorMatrix>(b.ptr, b.r, b.c);
    }

    BlockInfo& getBlockInfo(size_t index) {
        return blockInfos[index];
    }
};

}   // namespace gtsam
