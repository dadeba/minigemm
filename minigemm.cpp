#include <random>
#include <iostream>
#include <chrono>
#include "typetype.hpp"

#ifdef __linux__
//#include <gmpxx.h>
//#include <qd/qd_real.h>
//#include <quadmath.h>
#include "linux-perf-events.h"
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

#ifdef __linux__
  LinuxEvents<PERF_TYPE_HARDWARE> linux_events(std::vector<int>{
      PERF_COUNT_HW_CPU_CYCLES,
      PERF_COUNT_HW_INSTRUCTIONS,
  });
  std::vector<unsigned long long> results(2);
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

#ifdef __linux__    
    linux_events.start();
#endif
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
#ifdef __linux__
    linux_events.end(results);
    ccc += (double)results[0]; // cycles;
#endif
    
    time_after = std::chrono::steady_clock::now();
    double time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(time_after - time_before).count();
    elapsedtime += time_in_ns;
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
  double ops = 2*pow((double)n,3);
  printf("%5d %5d %5d %5.3e %e %e : %g GFlops\n", (int)m, (int)n, (int)k,
	 elapsedtime,  ccc, ccc/ops, ops/elapsedtime/1.0e9);

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
