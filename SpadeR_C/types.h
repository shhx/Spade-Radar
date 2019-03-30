#pragma once

#include <complex.h>

typedef struct data_mem {
	__int16* d_data;				// Read data
	size_t read_size;
	_Dcomplex* demodulated;	// DDC data
	_Dcomplex* filtered;		// LPF data
	size_t filtered_size;
	_Dcomplex* result;			// Matched filter data
	size_t result_length;
	_Dcomplex signal;
	fftw_complex *signal_FFT, *temp_result;	// FFT filter internal data
	_Dcomplex* temp_signal; //polyphasic filter data
}data_mem_t;

typedef struct config {
	size_t read_size;
	int decimate_factor;
	int fir_length;
	int matched_filt_deci_length;
	int nfft;
}config_t;

struct common_data {
	double* demod_table;
	fftw_complex *matched_filter_fft;
	fftw_plan p_forward, p_backward;
};