#include "FFT.h"

void FFT::fft(const std::vector<ComplexReal>& data, std::vector<ComplexReal>& output)
{
    int datasize = data.size();
    output.clear();
    output.resize(datasize);

    fftwf_complex *in;
    fftwf_complex *out;
    fftwf_plan plan;

    in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*datasize);
    out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*datasize);


    for (int i=0;i<datasize;i++)
    {
        in[i][0] = data[i].real();
        in[i][1] = data[i].imag();
    }

    plan = fftwf_plan_dft_1d(datasize,in,out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_execute(plan);

    for (int i=0;i<datasize;i++)
    {
        output[i] = {out[i][0], out[i][1]};
    }

    fftwf_destroy_plan(plan);
    fftwf_free(in);
    fftwf_free(out);

    return;
}

void FFT::ifft(const std::vector<ComplexReal>& data, std::vector<ComplexReal>& output)
{
    int datasize = data.size();

    output.clear();
    output.resize(datasize);

    fftwf_complex *in,*out;
    fftwf_plan plan;

    in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*datasize);
    out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*datasize);

    for (int i=0;i<datasize;i++)
    {
        in[i][0] = data[i].real();
        in[i][1] = data[i].imag();
    }

    plan = fftwf_plan_dft_1d(datasize, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftwf_execute(plan);

    for (int i=0;i<datasize;i++)
    {
        output[i] = {out[i][0], out[i][1]};
    }

    for (int i=0;i<datasize;i++)
    {
        output[i] /= datasize;
    }

    fftwf_destroy_plan(plan);
    fftwf_free(in);
    fftwf_free(out);

    return;
}

void FFT::autocorrelation(std::vector<std::vector<Real>>& data, std::vector<std::vector<Real>>& AC_vector)
{
    int dimension = data[0].size();
    int N = data.size();
    int N2 = N*2;

    // we have an autocorrelation for each dimension
    AC_vector.clear();
    AC_vector.resize(N, std::vector<Real>(dimension,0.0));

    for (int i=0;i<dimension;i++)
    {
        // perform forward fft
        std::vector<ComplexReal> forward_result;
        std::vector<ComplexReal> input(N2, {0,0});

        for (int j=0;j<N;j++)
        {
            input[j] = {data[j][i],0.0};
        }
        FFT::fft(input, forward_result);

        // obtain the Power spectral density (PSD)
        std::vector<ComplexReal> PSD(N2);
        for (int j=0;j<N2;j++)
        {
            Real square = std::pow(forward_result[j].real(),2.0) + std::pow(forward_result[j].imag(),2.0);
            PSD[j] = {square,0.0};
        }

        // perform inverse fft
        std::vector<ComplexReal> ifft;
        FFT::ifft(PSD, ifft);

        // copy the ifft result to autocorrelation vector
        for (int j=0;j<N;j++)
        {
            AC_vector[j][i] = ifft[j].real();
        }
    }
}