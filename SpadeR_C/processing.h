#pragma once
#include "utils.h"


void convolve(const double* signal, size_t signal_len, const double* h, size_t h_len, double* result);
void convolve_complex(const _Dcomplex* signal, size_t signal_len, const _Dcomplex* h, size_t h_len, _Dcomplex* result);
void filt_FIR(const _Dcomplex* signal, size_t signal_len, double* h, size_t h_len, _Dcomplex* result);
void filt_FIR_fast(const _Dcomplex* signal, size_t signal_len, const double* h, size_t h_len, _Dcomplex* result);
void filt_polyphasic_4(data_mem_t* data_mem, const _Dcomplex* signal, size_t signal_len, const double* h, size_t h_len, _Dcomplex* result);
void demodulate_IQ(const __int16* signal, size_t signal_len, const double* demod_table, _Dcomplex* result);
void generate_matched_filter(_Dcomplex* filter, size_t length, double pulse_width, double bandwidth, double FC, double FS, int real_filter, int decimate_factor);
void generate_demod_table(double* demod_table, double fc, double fs, size_t length);
double* decimate(const double* signal, size_t length_in, int factor);
_Dcomplex* decimate_complex(const _Dcomplex* signal, size_t length_in, int factor);
void convolve_overlap_add(const _Dcomplex* signal, size_t signal_len, const fftw_complex* h_FFT, size_t h_len, _Dcomplex* result, const fftw_plan p_forward, const fftw_plan p_backward, size_t nfft);
void convolve_fft(data_mem_t* data_mem, size_t signal_len, const fftw_complex* h_FFT, size_t h_len, _Dcomplex* result, const fftw_plan p_forward, const fftw_plan p_backward, size_t nfft);