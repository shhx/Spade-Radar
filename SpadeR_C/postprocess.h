#pragma once

typedef struct {
	_Dcomplex *data;
	size_t data_length;
	size_t intgrt_length;
	size_t index;
} matrix_t;

matrix_t* matrix_data_init(size_t data_length, size_t integration_length);
size_t matrix_data_add(matrix_t *m, const _Dcomplex *data, size_t data_length);
double* matrix_integrate(matrix_t *m);
size_t detect(double* signal, size_t signal_len, double* cfar, size_t* detected);
void cfar(double* signal, size_t signal_len, double* h_cfar, size_t cfar_length, double* result);