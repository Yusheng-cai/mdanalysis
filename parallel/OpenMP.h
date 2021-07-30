#pragma once
#include <omp.h>

namespace OpenMP 
{
    // Returns number of threads in the current team
    int get_num_threads();

    // Returns the thread number, with the current team, of the calling thread
    int get_thread_num();

    // obtain number of processors
    int get_num_procs();

    // obtain the maximum number of threads
    int get_max_threads();
} // namespace  
