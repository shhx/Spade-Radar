#include "init.h"

// Allocate all memory used by the CUDA kernels
void init_data_mem(data_mem_t* data_mem, size_t read_size, int decimate_factor, int fir_length, int matched_filt_deci_length, int nfft) {
	data_mem->read_size = read_size;
	data_mem->demodulated = (_Dcomplex*)malloc(read_size * sizeof(_Dcomplex));

	data_mem->filtered_size = read_size + fir_length - 1;
	data_mem->filtered = (_Dcomplex*)malloc(data_mem->filtered_size * sizeof(_Dcomplex));
	memset(data_mem->filtered, 0, data_mem->filtered_size * sizeof(_Dcomplex));

	data_mem->result_length = data_mem->filtered_size / decimate_factor + matched_filt_deci_length - 1;
	data_mem->result = (_Dcomplex*)malloc(data_mem->result_length * sizeof(_Dcomplex));
	memset(data_mem->result, 0, data_mem->result_length * sizeof(_Dcomplex));

	data_mem->signal_FFT = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nfft);
	data_mem->temp_result = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nfft);

	int D = 4;
	size_t length = (D*D + (D - 1) * 2 + read_size);
	data_mem->temp_signal = (_Dcomplex*)malloc(length * sizeof(_Dcomplex));
	memset(data_mem->temp_signal, 0, length * sizeof(_Dcomplex));
}