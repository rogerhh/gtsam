#include <stdint.h>
#include <stddef.h>
#include <assert.h>
#include <sched.h>
#include <unistd.h>
#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <time.h>
#ifndef BAREMETAL
#include <sys/mman.h>
#endif
#define NUM_CORE 4 // number of multithreading
#include "baremetal_tests/incremental_sphere2500_steps-2-200_period-25/incremental_dataset.h"

#include "cholesky.h"

pthread_barrier_t barrier_global;

#define INTERVAL 25

struct timespec got_node_start, got_node_end, AtA_start, AtA_end, chol_fac_end, chol_mem_end, chol_syrk_end, chol_end, merge_end;
struct timespec step_start, step_end;

float** node_workspaces = NULL;
pthread_mutex_t* node_locks;
int* node_num_children;
int* node_done_children;
int** node_children;

pthread_mutex_t queue_lock;
int node_ready_index = 0;
int node_ready_size = 0;
int* node_ready_queue;
int* node_affinity;
int num_active_nodes = 0;
// // record timing
// uint64_t* node_time; // TODO: AtA/cholesky breakdown 
// uint64_t* AtA_time;
// uint64_t* merge_time;
// uint64_t* chol_time;
// bool* chol_marked;
// uint64_t* init_time;
// uint64_t* backsolve_time;

typedef struct worker_args_t {
    int thread_id;
    int step;
} worker_args;

void* worker_cholesky(void* args_ptr) {
    worker_args* args = (worker_args*) args_ptr;
    int thread_id = args->thread_id;
    int step = args->step;
    int true_step = step + timestep_start;
    int nnodes = step_nnodes[step];

    bool* node_marked = step_node_marked[step];
    bool* node_fixed = step_node_fixed[step];

    int* node_parent = step_node_parent[step];
    int* node_height = step_node_height[step];
    int* node_width = step_node_width[step];
    float** node_data = step_node_data[step];
    int* node_num_blks = step_node_num_blks[step];
    int** node_A_blk_start = step_node_A_blk_start[step];
    int** node_B_blk_start = step_node_B_blk_start[step];
    int** node_blk_width = step_node_blk_width[step];

    int* node_num_factors = step_node_num_factors[step];
    int** node_factor_height = step_node_factor_height[step];
    int** node_factor_width = step_node_factor_width[step];
    float*** node_factor_data = step_node_factor_data[step];
    int** node_factor_num_blks = step_node_factor_num_blks[step];
    int*** node_factor_A_blk_start = step_node_factor_A_blk_start[step];
    int*** node_factor_B_blk_start = step_node_factor_B_blk_start[step];
    int*** node_factor_blk_width = step_node_factor_blk_width[step];


    clock_gettime(CLOCK_MONOTONIC, &got_node_start);

    double got_node_time = 0;
    double AtA_time = 0, chol_time = 0, merge_time = 0, chol_fac_time = 0, chol_mem_time = 0, chol_syrk_time = 0;    

    pthread_barrier_wait(&barrier_global);
    while(true) {
        
        bool done_flag = false;
        bool wait_flag = false;
        int node = -1;

        pthread_mutex_lock(&queue_lock);

        if(node_ready_index < node_ready_size) {
            node = node_ready_queue[node_ready_index];
            node_ready_index++;
        }
        else if(node_ready_size == num_active_nodes) {
            done_flag = true;
        }
        else if(node_ready_size > num_active_nodes) {
            printf("thread %d: Ready queue size greater than active nodes!\n", thread_id);
            exit(1);
        }
        else {
            wait_flag = true;
        }

        pthread_mutex_unlock(&queue_lock);

        if(done_flag) {
            break;
        }
        else if(wait_flag) {
            continue;
        }

        clock_gettime(CLOCK_MONOTONIC, &got_node_end);
        got_node_time += (double) (got_node_end.tv_nsec - got_node_start.tv_nsec) / 1000;
            
        // printf("thread %d node = %d\n", thread_id, node);

        bool marked = node_marked[node];
        bool fixed = node_fixed[node];

        if(!marked && !fixed) {
            printf("Node not marked or fixed!\n");
            exit(1);
        }

        // Technically we don't need to lock node because no two threads can
        // grab the same node. But lock it anyways

        // pthread_mutex_lock(&node_locks[node]);

        int parent = node_parent[node];
        int H_h = node_height[node];
        int H_w = node_width[node];
        float* H_data = node_data[node];

        int num_factors = node_num_factors[node];
        int* factor_height = node_factor_height[node];
        int* factor_width = node_factor_width[node];
        float** factor_data = node_factor_data[node];
        int* factor_num_blks = node_factor_num_blks[node];
        int** factor_A_blk_start = node_factor_A_blk_start[node];
        int** factor_B_blk_start = node_factor_B_blk_start[node];
        int** factor_blk_width = node_factor_blk_width[node];

        if(node_workspaces[node] == NULL) {
            // node_workspaces[node] = (float*) malloc(H_h * H_h * sizeof(float));
            node_workspaces[node] = (float*) my_malloc(H_h * H_h * sizeof(float));
            memset(node_workspaces[node], 0, H_h * H_h * sizeof(float));
        }

	
        float* ABC = node_workspaces[node];
        
        // Marked nodes
        // 1. AtA
        // 2. Cholesky (partial_factorization)
        // 3. Add LC to parent
        // 4. Copy [A B] back from workspace
        // Fixed nodes
        // 1. AtA
        // 2. LC = -LB LB^T
        // 3. Add LC to parent
        // 4. Don't copy [A B] back from workspace

        clock_gettime(CLOCK_MONOTONIC, &AtA_start);

        // 1. AtA
        for(int i = 0; i < num_factors; i++) {
            int h = factor_height[i];
            int w = factor_width[i];
            float* data = factor_data[i];

            int num_blks = factor_num_blks[i];
            int* A_blk_start = factor_A_blk_start[i];
            int* B_blk_start = factor_B_blk_start[i];
            int* blk_width = factor_blk_width[i];

            float* workspace = (float*) malloc(h * h * sizeof(float));
            memset(workspace, 0, h * h * sizeof(float));

            matmul(h, h, w, 
                   data, data, workspace,
                   h, h, h,
                   1, 1, 
                   true, false);


            sparse_matrix_add3_3(workspace, h, h, 
                                 ABC, H_h, H_h, 
                                 1, 
                                 num_blks, A_blk_start, B_blk_start, blk_width);

            free(workspace);
        }

        clock_gettime(CLOCK_MONOTONIC, &AtA_end);

        if(marked) {
            // Manually insert 1 in the diagonal
            float* ABC_col = ABC;
            for(int i = 0; i < H_w; i++) {
                *ABC_col = 1;
                ABC_col += H_h + 1;
            }

            // 2. Cholesky
            partial_factorization4(ABC, H_w, H_h);

            // 4. Copy [A B] back from workspace
            memcpy(H_data, ABC, H_w * H_h * sizeof(float));
        }
        else if(fixed) {
            // 2. LC = -LB LB^T
            int subdiag_h = H_h - H_w;
            float* LB = H_data + H_w;
            float* LC = ABC + H_w * (H_h + 1);

            matmul(subdiag_h, subdiag_h, H_w,
                   LB, LB, LC,
                   H_h, H_h, H_h,
                   -1, 1,
                   true, false);
        }
	// chol_marked[node] = marked ? true : false;

        clock_gettime(CLOCK_MONOTONIC, &chol_end);

        // 3. Add LC to parent
        if(parent != nnodes - 1 && parent != -1) {
            // printf("thread %d node = %d parent = %d\n", thread_id, node, parent);

            // lock parent
            pthread_mutex_lock(&node_locks[parent]);

            int subdiag_h = H_h - H_w;
            float* C = ABC + H_w * (H_h + 1);
            int next_H_h = node_height[parent];
            int next_H_w = node_width[parent];

            if(node_workspaces[parent] == 0) {
                // node_workspaces[parent] = (float*) malloc(next_H_h * next_H_h * sizeof(float));
                node_workspaces[parent] = (float*) my_malloc(next_H_h * next_H_h * sizeof(float));
                memset(node_workspaces[parent], 0, next_H_h * next_H_h * sizeof(float));
            }

            float* next_H_data = node_workspaces[parent];

            int num_blks = node_num_blks[node];
            int* A_blk_start = node_A_blk_start[node];
            int* B_blk_start = node_B_blk_start[node];
            int* blk_width = node_blk_width[node];

            sparse_matrix_add3_3(C, subdiag_h, H_h, 
                                 next_H_data, next_H_h, next_H_h, 
                                 1, 
                                 num_blks, A_blk_start, B_blk_start, blk_width);

            node_done_children[parent]++;
            if(node_done_children[parent] == node_num_children[parent]) {
                pthread_mutex_lock(&queue_lock);
                node_ready_queue[node_ready_size] = parent;
                node_affinity[node_ready_size] = thread_id;
                node_ready_size++;
                pthread_mutex_unlock(&queue_lock);
            }

            pthread_mutex_unlock(&node_locks[parent]);
        }

        clock_gettime(CLOCK_MONOTONIC, &merge_end);

        // 4. Free this nodes workspace
        // free(node_workspaces[node]);
        node_workspaces[node] = NULL;
        // pthread_mutex_unlock(&node_locks[node]);

        AtA_time += (double) ((long) (AtA_end.tv_nsec - AtA_start.tv_nsec)) / 1e6;
        chol_time += (double) ((long) (chol_end.tv_nsec - AtA_end.tv_nsec)) / 1e6;
        merge_time += (double) ((long) (merge_end.tv_nsec - chol_end.tv_nsec)) / 1e6;

        clock_gettime(CLOCK_MONOTONIC, &got_node_start);
    }

    printf("thread %d: AtA time = %f ms, chol time = %f ms, merge time = %f ms\n", thread_id, AtA_time, chol_time, merge_time);

    pthread_exit(NULL);
}

void* worker_backsolve(void* args_ptr) {
    worker_args* args = (worker_args*) args_ptr;
    int thread_id = args->thread_id;
    int step = args->step;
    int true_step = step + timestep_start;
    int nnodes = step_nnodes[step];

    bool* node_marked = step_node_marked[step];
    bool* node_fixed = step_node_fixed[step];

    int* node_parent = step_node_parent[step];
    int* node_height = step_node_height[step];
    int* node_width = step_node_width[step];
    int** node_ridx = step_node_ridx[step];
    float** node_data = step_node_data[step];
    float* x_data = step_x_data[step];

    pthread_barrier_wait(&barrier_global);
    while(true) {
        bool done_flag = false;
        bool wait_flag = false;
        int node = -1;
        pthread_mutex_lock(&queue_lock);

        if(node_ready_index < node_ready_size) {
            node = node_ready_queue[node_ready_index];
            node_ready_index++;
        }
        else if(node_ready_size == num_active_nodes) {
            done_flag = true;
        }
        else if(node_ready_size > num_active_nodes) {
            printf("thread %d: Ready queue size greater than active nodes!\n", thread_id);
            exit(1);
        }
        else {
            wait_flag = true;
        }

        pthread_mutex_unlock(&queue_lock);

        if(done_flag) {
            break;
        }
        else if(wait_flag) {
            continue;
        }
            
        //printf("thread %d solve node = %d\n", thread_id, node);
        int width = node_width[node];
        int height = node_height[node];
        int* ridx = node_ridx[node];

        partial_backsolve(node_data[node], width, height, height, ridx, x_data);

        pthread_mutex_lock(&queue_lock);
        for(int i = 0; i < node_num_children[node]; i++) {
            int child = node_children[node][i];
            node_ready_queue[node_ready_size] = child;
            node_ready_size++;

        }
        pthread_mutex_unlock(&queue_lock);

    }

    pthread_exit(NULL);
}

void *print_message(void *ptr){
    int cpu_id = sched_getcpu();
   // char *msg;
   // msg = (char *) ptr;
    printf("print msg - cpu_id: %d \n", cpu_id);
   // printf("%s \n", msg);
   return NULL;
}

int main(int argc, char** argv) {
    int num_threads = atoi(argv[1]);

    int cpu_id;
    cpu_id = sched_getcpu();
    printf("main thread cpuid: %d \n", cpu_id);

    cpu_set_t cpuset[NUM_CORE];
    pthread_t thread[NUM_CORE];
    pthread_attr_t attr[NUM_CORE];
    for(int i = 0; i < NUM_CORE; i++)
      pthread_attr_init(&attr[i]);

    printf("create threading \n");
    for(int i = 0; i < num_threads; i++){
	  CPU_ZERO(&cpuset[i]);
	  CPU_SET(i, &cpuset[i]);
	  pthread_attr_setaffinity_np(&attr[i], sizeof(cpu_set_t), &cpuset[i]);
	  pthread_create(&thread[i], &attr[i], print_message, NULL);
    }

    for(int i = 0; i < num_threads; i++){
      pthread_join(thread[i], NULL);
    }
    printf("thread joined after message printing\n");

    pthread_barrier_init(&barrier_global, NULL, num_threads);
    for(int step = 0; step < num_timesteps; step++) {
        clock_gettime(CLOCK_MONOTONIC, &step_start);

        int true_step = (step+1)*INTERVAL;//step + timestep_start;
        printf("true step = %d\n", true_step);

        int nnodes = step_nnodes[step];
        // node_time = (uint64_t*) malloc(nnodes*sizeof(uint64_t));
        // memset(node_time, 0, nnodes*sizeof(uint64_t));
        // init_time = (uint64_t*) malloc(nnodes*sizeof(uint64_t));
        // AtA_time = (uint64_t*) malloc(nnodes*sizeof(uint64_t));
        // merge_time = (uint64_t*) malloc(nnodes*sizeof(uint64_t));
        // chol_time = (uint64_t*) malloc(nnodes*sizeof(uint64_t));
        // chol_marked = (bool*) malloc(nnodes*sizeof(bool));
        // backsolve_time = (uint64_t*) malloc(nnodes*sizeof(uint64_t));
        // memset(backsolve_time, 0, nnodes*sizeof(uint64_t));

        bool* node_marked = step_node_marked[step];
        bool* node_fixed = step_node_fixed[step];

        int* node_parent = step_node_parent[step];
        int* node_height = step_node_height[step];
        int* node_width = step_node_width[step];
        float** node_data = step_node_data[step];
        int* node_num_blks = step_node_num_blks[step];
        int** node_A_blk_start = step_node_A_blk_start[step];
        int** node_B_blk_start = step_node_B_blk_start[step];
        int** node_blk_width = step_node_blk_width[step];

        int* node_num_factors = step_node_num_factors[step];
        int** node_factor_height = step_node_factor_height[step];
        int** node_factor_width = step_node_factor_width[step];
        float*** node_factor_data = step_node_factor_data[step];
        int** node_factor_num_blks = step_node_factor_num_blks[step];
        int*** node_factor_A_blk_start = step_node_factor_A_blk_start[step];
        int*** node_factor_B_blk_start = step_node_factor_B_blk_start[step];
        int*** node_factor_blk_width = step_node_factor_blk_width[step];
        
        // node_num_children[node] gives the number of active children
        // in this node. When node_done_childre[node] == node_num_children[node]
        // then node is available to work on
        // Each node_locks[node] controls node_done_children[node] and node_workspaces[node]
        // queue_lock controls node_ready_queue
        pthread_mutex_init(&queue_lock, NULL);

        node_locks = (pthread_mutex_t*) malloc(nnodes * sizeof(pthread_mutex_t));
        node_num_children = (int*) malloc(nnodes * sizeof(int));
        node_done_children = (int*) malloc(nnodes * sizeof(int));
        node_children = (int**) malloc(nnodes * sizeof(int*));
        node_ready_queue = (int*) malloc(nnodes * sizeof(int));
        node_affinity = (int*) malloc(nnodes * sizeof(int));
        memset(node_num_children, 0, nnodes * sizeof(int));
        memset(node_done_children, 0, nnodes * sizeof(int));
        memset(node_children, 0, nnodes * sizeof(int*));

        num_active_nodes = 0;
        node_ready_size = 0;
        node_ready_index = 0;

        for(int node = 0; node < nnodes - 1; node++) {
            bool marked = node_marked[node];
            bool fixed = node_fixed[node];

            pthread_mutex_init(&node_locks[node], NULL);

            if(!marked && !fixed) { continue; }

            num_active_nodes++;

            if(node_num_children[node] == 0) {
                node_ready_queue[node_ready_size] = node;
                node_affinity[node_ready_size] = 0;
                node_ready_size++;
            }

            int parent = node_parent[node];
            node_num_children[parent]++;
        }
/*
        printf("active node size: %d\n", num_active_nodes);
        printf("num children: \n");
        for(int node = 0; node < nnodes - 1; node++) {
            printf("%d ", node_num_children[node]);
        }
        printf("\n");
*/
        // This is used as workspaces for each node's ABC matrix
        // We will allocate a node's workspace when its first child
        // is done
        node_workspaces = (float**) malloc(nnodes * sizeof(float*));
        memset(node_workspaces, 0, nnodes * sizeof(float*));

        pthread_t threads[NUM_CORE];
        worker_args args[NUM_CORE];
        for(int thread = 0; thread < num_threads; thread++) {
            args[thread].thread_id = thread;
            args[thread].step = step;
            pthread_create(&threads[thread], &attr[thread], worker_cholesky, &args[thread]);
        }

        for(int thread = 0; thread < num_threads; thread++) {
            pthread_join(threads[thread], NULL);
        }

        // Reset the ready queue for backsolve
        node_ready_index = 0;
        node_ready_size = 0;

        // Backsolve. First set up a list of children for each node
        for(int node = nnodes - 2; node >= 0; node--) {
            bool marked = node_marked[node];
            bool fixed = node_fixed[node];

            if(!marked && !fixed) { continue; }

            node_children[node] = (int*) malloc(node_num_children[node] * sizeof(int));

            // Reset this so we can have some sort of indexing
            node_num_children[node] = 0;

            int parent = node_parent[node];
            if(parent == nnodes - 1) {
                node_ready_queue[node_ready_size] = node;
                node_ready_size++;
            }
            else {
                node_children[parent][node_num_children[parent]] = node;
                node_num_children[parent]++;
            }
        }
        for(int thread = 0; thread < num_threads; thread++) {
            args[thread].thread_id = thread;
            args[thread].step = step;
            pthread_create(&threads[thread], &attr[thread], worker_backsolve, &args[thread]);
        }

        for(int thread = 0; thread < num_threads; thread++) {
            pthread_join(threads[thread], NULL);
        }

        for(int node = 0; node < nnodes - 1; node++) {
            pthread_mutex_destroy(&node_locks[node]);
        }
        pthread_mutex_destroy(&queue_lock);

        for(int node = 0; node < nnodes - 1; node++) {
            free(node_children[node]);
        }

        free(node_locks);
        free(node_num_children);
        free(node_done_children);
        free(node_children);
        free(node_ready_queue);
        free(node_workspaces);
	my_free_all(NULL);

        node_workspaces = NULL;
        node_ready_queue = NULL;
        node_done_children = NULL;
        node_num_children = NULL;
        node_children = NULL;
        node_locks = NULL;
  
        int nnode = step_nnodes[step];

        // free(merge_time);
        // free(chol_time);
        // free(chol_marked);
        // free(init_time);
        // free(node_time);
        // free(AtA_time);
        // free(backsolve_time);
        // merge_time = NULL;
        // chol_time = NULL;
        // init_time = NULL;
	// chol_marked = NULL;
        // node_time = NULL;
        // AtA_time = NULL;
        // backsolve_time = NULL;

	clock_gettime(CLOCK_MONOTONIC, &step_end);
	double cpu_time_used = (step_end.tv_sec - step_start.tv_sec) * 1000
            + ((long) (step_end.tv_nsec - step_start.tv_nsec)) / 1000000;
	printf("step %d time: %f ms\n", step, cpu_time_used);

    }
    pthread_barrier_destroy(&barrier_global);
    printf("end of test\n");
}
