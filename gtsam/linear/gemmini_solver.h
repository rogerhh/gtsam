#pragma once

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdint.h>

typedef struct lls_solver_t {
  // Global variables that do not change step to step go here
  int num_threads;
} lls_solver;

typedef struct lls_solver_args_t {
  // Step-dependent variables go here
  bool no_values;
  int step;
  int nnodes;

  bool* node_marked;
  bool* node_fixed;

  int* node_parent;
  int* node_height;
  int* node_width;
  float** node_data;
  int* node_num_blks;
  int** node_A_blk_start;
  int** node_B_blk_start;
  int** node_blk_width;
  int** node_ridx;
  float* x_data;

  int* node_num_factors;
  int** node_factor_height;
  int** node_factor_width;
  float*** node_factor_data;
  int** node_factor_num_blks;
  int*** node_factor_A_blk_start;
  int*** node_factor_B_blk_start;
  int*** node_factor_blk_width;

} lls_solver_args;

void lls_solver_init(lls_solver* solver);
void lls_solver_destroy(lls_solver* solver);

void lls_solver_solve(lls_solver* solver, lls_solver_args* solver_args);



