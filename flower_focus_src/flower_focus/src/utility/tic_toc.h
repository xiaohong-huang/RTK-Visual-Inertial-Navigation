

#pragma once

#include <ctime>
#include <cstdlib>
#include <chrono>

#include <omp.h>


class TicToc {
  public:
    TicToc() {
        tic();
    }

    void tic() {
        start = omp_get_wtime();
    }

    double toc() {
        end = omp_get_wtime();
        return (end - start) * 1000;
    }

    double toc_s() {
        end = omp_get_wtime();
        return end - start;
    }

  private:
    double start, end;
};
