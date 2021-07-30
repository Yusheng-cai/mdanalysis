#include <mpi.h>

namespace MPI_environment {
    // initialize MPI
    void init(int argc, char** argv);

    // check if MPI is initialized
    bool is_initialized();

    // finalize MPI
    void finalize();
}