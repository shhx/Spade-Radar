#pragma once

#include "stdafx.h"

void init_data_mem(data_mem_t* data_mem, size_t read_size, int decimate_factor, int fir_length, int matched_filt_deci_length, int nfft);