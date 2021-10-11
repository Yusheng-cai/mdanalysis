#include "fftw3.h"
#include "tools/CommonTypes.h"

#include <vector>
#include <array>
#include <string>
#include <complex>
#include "parallel/OpenMP.h"

namespace FFT
{
    using Real = CommonTypes::Real;
    using ComplexReal = std::complex<Real>;

    void fft(const std::vector<ComplexReal>& data, std::vector<ComplexReal>& output);
    void ifft(const std::vector<ComplexReal>& data, std::vector<ComplexReal>& output);

    void autocorrelation(std::vector<std::vector<Real>>& data, std::vector<std::vector<Real>>& AC_vector, bool biased);
}
