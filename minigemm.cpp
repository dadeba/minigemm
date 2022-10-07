#include <random>
#include <iostream>
#include <chrono>
#include "typetype.hpp"
#include "m1cycles.h"

#ifdef __linux__
#include "linux-perf-events.h"
#endif

#ifdef HIGHPREC
#include <gmpxx.h>
#include <qd/qd_real.h>
#endif
  
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

  performance_counters agg_min{1e300};
  performance_counters agg_avg{0.0};
  
  double elapsedtime = 0.0;
  const int LOOP = 32;
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

    performance_counters start = get_counters();
    for (int j = 0; j < n; j++) {
      for (int l = 0; l < k; l++) {
	REAL temp = b[j*lda+l];
	for (int i = 0; i < m; i++) {
	  c[j*lda+i] += temp * a[l*lda+i];
	}
      }
    }
    performance_counters end = get_counters(true);
    performance_counters diff = end - start;

    agg_min = agg_min.min(diff);
    agg_avg += diff;
    
    time_after = std::chrono::steady_clock::now();
    double time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(time_after - time_before).count();
    elapsedtime += time_in_ns;
  }
  elapsedtime = elapsedtime * 1.0e-9 / (double)LOOP;

  agg_avg /= (double)LOOP;
  agg_min /= 1;

  //  printf("    m     n     k     sec    cycles\n");
  double ops = 2*pow((double)n,3);
  printf("%5d : %8.3e %3.2e(%3.2e) %3.2e(%3.2e) : %3.2e : %8.4g GFlops\n", (int)n,
	 elapsedtime,
	 agg_avg.cycles, agg_min.cycles,
	 agg_avg.instructions, agg_min.instructions,
	 agg_avg.cycles/ops, ops/elapsedtime/1.0e9);

  delete[] c;
  delete[] b;
  delete[] a;
}

int main(int argc, char *argv[])
{
  setup_performance_counters();
  
  const int nn = 199;
  minigemm<float>(nn);
  minigemm<double>(nn);
  minigemm<_Float16>(nn);

#ifdef HIGHPREC
  minigemm<dd_real>(nn);
  minigemm<qd_real>(nn);
  minigemm<mpf_class>(nn);
  minigemm<_Float128>(nn);
#endif
}
