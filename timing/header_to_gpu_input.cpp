#include "baremetal_tests/incremental_sphere2500_steps-2-200_period-25/incremental_dataset.h"

#include <set>
#include <map>

using namespace std;

int main() {

    // Copy all the factors out of the dataset
    // and remap all the factor indices

    for(int step = 0; step < num_timesteps; step++) {
        int true_step = step + timestep_start;

        int nnodes = step_nnodes[step];

        bool* node_marked = step_node_marked[step];
        bool* node_fixed = step_node_fixed[step];

        int** node_ridx = step_node_ridx[step];

        int* node_num_factors = step_node_num_factors[step];
        int** node_factor_height = step_node_factor_height[step];
        int** node_factor_width = step_node_factor_width[step];
        int*** node_factor_ridx = step_node_factor_ridx[step];

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
                        h_csrColIndA.push_back(node_ridx[node][ridx[ih]]);
                        if(ih == j) { h_csrValA.push_back(1.0f); }
                        else { h_csrValA.push_back(0.0f); }
                    }
                }
            }
        }
    }
}
