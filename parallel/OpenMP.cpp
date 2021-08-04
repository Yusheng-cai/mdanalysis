#include "OpenMP.h"

int OpenMP::get_num_threads()
{
    return omp_get_num_threads();
}

int OpenMP::get_thread_num()
{
    return omp_get_thread_num();
}

int OpenMP::get_num_procs()
{
    return omp_get_num_procs();
}

int OpenMP::get_max_threads()
{
    return omp_get_max_threads();
}

bool OpenMP::in_parallel()
{
    return omp_in_parallel();
}