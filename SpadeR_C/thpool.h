#pragma once

#include <windows.h>
#include <stdio.h>                                                              
#include <stdlib.h>                                                             
#include <errno.h>                                                              
#include <string.h> 
#include "init.h"

typedef struct tpool_work {
	void(*routine)();
	void *arg1;
	size_t arg2;
	struct tpool_work *next;
} tpool_work_t;

typedef struct tpool {
	int num_threads;
	int max_queue_size;

	int do_not_block_when_full;
	HANDLE *threads;
	int cur_queue_size;
	tpool_work_t *queue_head;
	tpool_work_t *queue_tail;
	HANDLE queue_lock;
	HANDLE queue_not_empty;
	HANDLE queue_not_full;
	HANDLE queue_empty;
	int queue_closed;
	int shutdown;
	config_t *config;
	struct common_data *common_data;
	data_mem_t *data_mem_th;
	DWORD *thIDs;
} *tpool_t;

void tpool_init(tpool_t *tpoolp, int num_worker_threads, int max_queue_size, int do_not_block_when_full, config_t* config, struct common_data *common_data);
int tpool_add_work(tpool_t tpool, void(*routine)(), void *arg1, size_t arg2);
int tpool_destroy(tpool_t tpoolp, int finish);

