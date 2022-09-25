#include <random>
#include <iostream>
#include <chrono>
#include "typetype.hpp"

#ifdef __linux__
//#include <gmpxx.h>
//#include <qd/qd_real.h>
//#include <quadmath.h>
#define RDTSC(X)  asm volatile ("rdtsc; shlq $32,%%rdx; orq %%rdx,%%rax" : "=a" (X) :: "%rdx")
#endif

#ifdef __APPLE__
#include "m1cycles.h"
void pretty_print(std::pair<performance_counters, performance_counters> result);
#define RDTSC(X)  do {} while(0);
#endif

#define MFLOPS 1e-6
#define NANOSECOND 1e-9
  
template <typename REAL>
void minigemm(const int nn)
{
  std::cout << "Test for " << TypeName<REAL>() << "\n";

  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<> urdist(0.0, 1.0);

  const int lda = nn;
  const int ldb = nn;
  const int ldc = nn;

  const int n = lda;
  const int k = lda;
  const int m = lda;
  
  REAL *a = new REAL [lda * n];
  REAL *b = new REAL [lda * n];
  REAL *c = new REAL [lda * n];

#ifdef __APPLE__  
  performance_counters agg_min{1e300};
  performance_counters agg_avg{0.0};
#endif

  double ccc = 0.0; // average cylces
  double elapsedtime = 0.0;
  const int LOOP = 1024;
  for (int jj = 0; jj < LOOP; jj++) {
    for(int j = 0; j < n; j++) {
      for(int i = 0; i < n; i++) {
	c[j*lda + i] = 0.0;
	a[j*lda + i] = urdist(mt);
	b[j*lda + i] = urdist(mt);
      }
    }

    std::chrono::steady_clock::time_point time_before;
    std::chrono::steady_clock::time_point time_after;

    time_before = std::chrono::steady_clock::now();   

    unsigned long long t1 = 0;
    unsigned long long t2 = 0;
    
    RDTSC(t1);
#ifdef __APPLE__      
    performance_counters start = get_counters();
#endif
    for (int j = 0; j < n; j++) {
      for (int l = 0; l < k; l++) {
	REAL temp = b[j*lda+l];
	for (int i = 0; i < m; i++) {
	  c[j*lda+i] += temp * a[l*lda+i];
	}
      }
    }
#ifdef __APPLE__
    performance_counters end = get_counters();
    performance_counters diff = end - start;
    agg_min = agg_min.min(diff);
    agg_avg += diff;
#endif
    RDTSC(t2);
    time_after = std::chrono::steady_clock::now();
    double time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(time_after - time_before).count();

    elapsedtime += time_in_ns;
    ccc += (double)(t2 - t1);
  }
  elapsedtime = elapsedtime * NANOSECOND / (double)LOOP;
  ccc = ccc/(double)LOOP;

#ifdef __APPLE__
  agg_avg /= (double)LOOP;
  agg_min /= 1;
  ccc = agg_avg.cycles;
  //  double ccc_min = agg_min.cycles;
#endif
  
  //  printf("    m     n     k     sec    cycles\n");
  //  printf("%5d %5d %5d %5.3f %e(%5.3f) %e\n", (int)m, (int)n, (int)k, elapsedtime, ccc, ccc_min, ccc/(2.0*n*n*n));
  printf("%5d %5d %5d %5.3f %e\n", (int)m, (int)n, (int)k, ccc, ccc/(2.0*n*n*n));
  delete[] c;
  delete[] b;
  delete[] a;
}

int main(int argc, char *argv[])
{
#ifdef __APPLE__
  setup_performance_counters();
#endif
  
  const int nn = 64;
  minigemm<float>(nn);
  minigemm<double>(nn);
  minigemm<_Float16>(nn);

#ifdef __linux__
  //  minigemm<dd_real>(nn);
  //  minigemm<qd_real>(nn);
  //  minigemm<mpf_class>(nn);
  minigemm<_Float128>(nn);
#endif
}

#ifdef __APPLE__
void pretty_print(std::pair<performance_counters, performance_counters> result)
{
  printf(" %32s ", "");
  printf(" %8.2f instructions/float (+/- %3.1f %%) ",
	 result.first.instructions,
         (result.second.instructions - result.first.instructions) * 100.0 /
	 result.first.instructions);

  printf("\n");
  printf(" %32s ", "");
  printf(" %8.2f cycles/float (+/- %3.1f %%) ", result.first.cycles,
         (result.second.cycles - result.first.cycles) * 100.0 /
             result.first.cycles);

  printf("\n");
  printf(" %32s ", "");
  printf(" %8.2f instructions/cycle ",
         result.first.instructions / result.first.cycles);
  printf("\n");
  printf(" %32s ", "");
  printf(" %8.2f branches/float (+/- %3.1f %%) ", result.first.branches,
         (result.second.branches - result.first.branches) * 100.0 /
             result.first.branches);

  printf("\n");
  printf(" %32s ", "");
  printf(" %8.4f mis. branches/float ", result.second.missed_branches);
  printf("\n");
}
#endif
