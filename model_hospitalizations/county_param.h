#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <math.h>
#include <time.h>
#include <unistd.h>
#include <sstream>
#include <cstring>
#include <string>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <algorithm>  // std::sort

struct CountyParams{
	int n_pop;
	int *pop_N;
	double *dist_vec;
	double *census_area;
	int *E_pops;
	int *dist_ID;
};


struct OutputParams{
	int *out_pops_seed; //not_used (for testing)
	int *out_pops_time; //not_used (for testing)
	int *out_pops_pop; //not_used (for testing)
	int *out_total_hosp;
    int *out_I_hosp;
    int *out_I_icu1;
    int *out_I_icu2;
    int *out_n_death;
	int n_pop;
};
