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

void FFT::autocorrelation(std::vector<std::vector<Real>>& data, std::vector<std::vector<Real>>& AC_vector, bool biased)
{
    int dimension = data[0].size();

    // we have an autocorrelation for each dimension
    AC_vector.clear();
    AC_vector.resize(dimension);

    int N = data.size();
    int N2 = N*2;

    for (int i=0;i<dimension;i++)
    {
        // to calculate autocorrelation, we must split the data up
        std::vector<Real> tempdata(N2,0.0);

        // resize the ith dimension AC to be data size large
        AC_vector[i].resize(N,0.0);
        
        for (int j=0;j<N;j++)
        {
            tempdata[j] = data[j][i];
        }

        std::vector<ComplexReal> fft;
        std::vector<ComplexReal> input(N2);

        for (int j=0;j<N2;j++)
        {
            ComplexReal number(tempdata[j],0);
            input[j] = number;
        }

        FFT::fft(input, fft);

        std::vector<ComplexReal> squared(N2);
        for (int j=0;j<N2;j++)
        {
            Real square = std::pow(fft[j].real(),2.0) + std::pow(fft[j].imag(),2.0);
            ComplexReal number(square,0.0);
            squared[j] = number;
        }

        std::vector<ComplexReal> ifft;
        FFT::ifft(squared, ifft);

        // divide by N or (N-m)
        if ( ! biased)
        {
            for (int j=0;j<N;j++)
            {
                AC_vector[i][j] = ifft[j].real()/N;
            }
        }
        else
        {
            for (int j=0;j<N;j++)
            {
                AC_vector[i][j] = ifft[j].real()/(N-j);
                std::cout << "AC vector " << j << " = " << AC_vector[i][j] << "\n";
            }
        }
    }
}