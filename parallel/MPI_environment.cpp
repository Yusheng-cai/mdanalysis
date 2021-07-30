#include "MPI_environment.h"

void MPI_environment::init(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
}

bool MPI_environment::is_initialized()
{
    int is_initialized = 0;
    MPI_Initialized(&is_initialized);

    if (is_initialized == 0){return false;}
    else {return true;}
}

void MPI_environment::finalize()
{
    MPI_Finalize();
}