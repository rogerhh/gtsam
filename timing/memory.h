#pragma once

#include <stdint.h>
#include <stdio.h>

#define MEMORY_SIZE 1638400 // bytes


// Define a char array to represent the memory space
static double _memory[MEMORY_SIZE / 8];
static int stack_ptr = 0;

// Function to allocate memory
void* my_malloc(unsigned int size) {
  
  
  // Check memory bound
  if(stack_ptr + size >= MEMORY_SIZE) {
    exit(1);
  }

  int remainder = size & 0b111;
  int size8 = remainder == 0? size / 8 : size / 8 + 1;
  double* ret = _memory + stack_ptr;
  stack_ptr += size8;

  return (void*) ret;
}

// Function to reallocate memory
void* my_realloc(void* ptr, int size);

// Function to free all allocated memory after ptr
void my_free_all(void* ptr) {
  stack_ptr = 0;
}

// Function to free all allocated memory after ptr
void my_free_after(void* ptr) {
  stack_ptr = (double*) ptr - _memory;
}

void my_memcpy(void* dest, void* src, unsigned int len) {
  char* d = (char*) dest;
  char* s = (char*) src;
  for(int i = 0; i < len; i++) {
    *d = *s;
    d++;
    s++;
  }
}

void my_memset(void* dest, char a, unsigned int len) {
  char* d = (char*) dest;
  for(int i = 0; i < len; i++) {
    *d = a;
    d++;
  }
}
