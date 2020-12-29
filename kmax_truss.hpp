#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <fcntl.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <time.h>
#include <omp.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <parallel/algorithm>

using namespace std;

int64_t get_filesize(char const *filename);

uint32_t* import_and_trans(char const *dt_name, size_t *ntmp);

void trans_txt(char *rawdata, uint32_t *endpoints, size_t dt_size, size_t n_edge);

void par_sort(uint64_t *ptrz, size_t n_edge);

void cal_ptfreq(uint32_t *endpoints, uint32_t *ntr_mark, uint32_t *freq_count, uint32_t n_node, size_t n_edge, uint32_t n_trust);

int64_t edge_prune(uint32_t *endpoints, uint32_t *ntr_mark, uint32_t *freq_count, uint32_t n_node, size_t n_edge, uint32_t n_trust);

int64_t* mark_rowhead(uint32_t *endpoints, uint32_t n_node, int64_t n_edge);

uint32_t count_tris(uint32_t *val_mx, uint32_t *pts_bak, int64_t n_ptbk, uint32_t n_node, uint32_t n_try);

int64_t ktruss_chk(uint32_t *endpoints, uint32_t n_node, size_t n_edge, uint32_t n_try);

uint32_t data_dim_zip(uint32_t *endpoints, uint32_t* freq_count, uint32_t n_node, size_t n_edge);

// uint32_t data_deep_zip(uint32_t *endpoints, uint32_t n_node, size_t n_edge);

// bool cmp(uint64_t a, uint64_t b);