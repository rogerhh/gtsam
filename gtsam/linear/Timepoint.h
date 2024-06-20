#pragma once

#include <chrono>

// #define NO_TIMING

#ifdef RISCV
inline uint64_t read_cycles() {
  uint64_t cycles;
  asm volatile ("rdcycle %0" : "=r" (cycles));
  return cycles;
}
#endif

struct Timepoint {
#ifndef NO_TIMING
#ifdef RISCV
  uint64_t t;
  Timepoint() : t(read_cycles()) {}
#else
  std::chrono::time_point<std::chrono::high_resolution_clock> t;
  Timepoint() : t(std::chrono::high_resolution_clock::now()) {}
#endif
#endif
};

inline int64_t operator-(const Timepoint& t1, const Timepoint& t2) {
#ifndef NO_TIMING
#ifdef RISCV
  return t1.t - t2.t;
#else
  return std::chrono::duration_cast<std::chrono::nanoseconds>(t1.t - t2.t).count();
#endif
#else
  return 0;
#endif
}

