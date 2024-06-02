#include "baremetal_tests/incremental_sphere2500_steps-2-200_period-25/incremental_dataset.h"

#include <iostream>
#include <vector>
#include <set>
#include <map>

using namespace std;

int main() {

    // Copy all the factors out of the dataset
    // and remap all the factor indices

    int num_total_factors = 0;

    cout << num_timesteps << endl;

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

        // Host data
        vector<int> h_csrRowPtrA = {0};
        vector<int> h_csrColIndA;
        vector<int> h_csrRowPtrAT = {0};
        vector<int> h_csrColIndAT;
        vector<float> h_b;
        vector<float> h_x;
        vector<int> h_csrRowPtrD;
        vector<int> h_csrColIndD;
        vector<int> h_csrRowPtrH;
        vector<int> h_csrColIndH;

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

            num_total_factors += num_factors;

            for(int i = 0; i < num_factors; i++) {
                int height = factor_height[i];
                int width = factor_width[i];
                int* ridx = factor_ridx[i];

                for(int ih = 0; ih < height - 1; ih++) {
                    ridx_set.insert(node_ridx[node][ridx[ih]]);
                }
            }
        }

        int count = 0;
        for(int ridx : ridx_set) {
            remapped_ridx[ridx] = count++;
        }

        cout << num_total_factors << endl;

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

                cout << height << " " << width << endl;

                for(int ih = 0; ih < height - 1; ih++) {
                    cout << remapped_ridx.at(node_ridx[node][ridx[ih]]) << " ";
                }
                cout << endl;
            }
        }
    }
}
