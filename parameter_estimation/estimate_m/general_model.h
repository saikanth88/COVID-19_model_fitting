#ifndef GENERAL_MODEL_H
#define GENERAL_MODEL_H

struct GeneralParamStruct
{
    int *input_realz_seeds;
    double *input_census_area;
    double *input_dist_vec;

    double *input_r0;           // vector from time windows
    double *input_dist_param;   // vector from time windows
    double *input_m;            // vector from time windows
    double *input_imm_frac;     // vector from time windows
    int *input_window_length;   // vector from time windows
    int total_windows;          // total number of time windows

    int n_realz;
    int n_pop;
    int t_max;
    int tau;            // remove?
    int n_equations;    // remove?
    int trans_type;
    double *beta;
    double m;
    double imm_frac;
    double stoch_sd;
    double dd_trans_monod_k;

    float *pop_N;
    float *census_area;
    float **dist_mat;
    float **prob_move;
};


#endif // GENERAL_MODEL_H
