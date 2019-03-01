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
