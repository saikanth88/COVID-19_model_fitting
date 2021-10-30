//Setting parameters for the covid19_model()

#include "covid19_model.cpp"

using namespace std;

OutputParams covid19_driver(int t_max, int n_realz, double *beta_vec, double *mgrt_vec, CountyParams par) {
    double stoch_sd = 0.05; //OK
    int trans_type = 1; //OK
    double dd_trans_monod_k = 500; //OK
    double frac_beta_asym = 0.55; //OK
    double frac_beta_hosp = 0.05; //OK
    double delta = 1/3.0; //OK
    double recov_a = 1/6.0; //OK
    double recov_p = 1/2.0; //OK
    double recov_s = 1/6.0; //OK
    double recov_home = 1/3.0; //OK
    double recov_icu1 = 1/8.0; //OK
    double recov_icu2 = 1/4.0; //OK
    double asym_rate = 0.28; //OK
    double sym_to_icu_rate = 0.015; //OK

    int n_pop = par.n_pop; //OK

    int nrow = n_pop * n_realz * t_max;

    vector<int> out_pops_seed(nrow);
    vector<int> out_pops_pop(nrow);
    vector<int> out_pops_time(nrow);
    vector<int> out_pops_S_pop(nrow);
    vector<int> out_pops_E_pop(nrow);
    vector<int> out_pops_I_asym_pop(nrow);
    vector<int> out_pops_I_presym_pop(nrow);
    vector<int> out_pops_I_sym_pop(nrow);
    vector<int> out_pops_I_home_pop(nrow);
    vector<int> out_pops_I_hosp_pop(nrow);
    vector<int> out_pops_I_icu1_pop(nrow);
    vector<int> out_pops_I_icu2_pop(nrow);
    vector<int> out_pops_R_pop(nrow);
    vector<int> out_pops_D_pop(nrow);
    vector<int> out_events_pos(nrow);
    vector<int> out_events_sym(nrow);
    vector<int> out_events_total_hosp(nrow);
    vector<int> out_events_total_icu(nrow);
    vector<int> out_events_n_death(nrow);

    vector<int> input_realz_seeds(n_realz);
    std::iota (input_realz_seeds.begin(), input_realz_seeds.end(), 1);

    vector<double> input_r0(t_max, 1.05); //not_used

    vector<double> input_dist_param(t_max, 25.0);
//    std::copy ( dist_vec, dist_vec + t_max, input_dist_param.begin());
    vector<double> input_m(t_max, 0.0);
    std::copy ( mgrt_vec, mgrt_vec + t_max, input_m.begin());

    vector<double> input_imm_frac(t_max, 0.0);

    if (t_max > 54){
    	std::fill (input_dist_param.begin()+54, input_dist_param.end(), 75.0);
//    	std::fill (input_m.begin()+54, input_m.end(), 0.03);
//    	std::fill (input_imm_frac.begin()+54, input_imm_frac.end(), 0.01);
   }

    vector<double> input_hosp_rate(t_max, 0.2); //OK
    vector<double> input_icu_rate(t_max, 0.28); //OK
    vector<double> input_death_rate(t_max, 0.6); //OK
    vector<double> input_recov_hosp(t_max, 1/7.0); //OK
    vector<int> input_window_length(t_max, 1);; //window length is 1 (daily calc)

    vector<double> input_census_area(n_pop, 0.0);
    std::copy ( par.census_area, par.census_area + n_pop, input_census_area.begin());

    vector<double> input_dist_mat(n_pop*n_pop, 0.0);
    std::copy ( par.dist_vec, par.dist_vec + n_pop*n_pop, input_dist_mat.begin());

    vector<int> input_N_pops(n_pop, 0);
    std::copy ( par.pop_N, par.pop_N + n_pop, input_N_pops.begin());

    vector<int> input_S_pops(n_pop, 0);

    vector<int> input_E_pops(n_pop, 0);
    std::copy ( par.E_pops, par.E_pops + n_pop, input_E_pops.begin());

    for(int k = 0; k < n_pop; k++)
    	input_S_pops[k] =  input_N_pops[k] -  input_E_pops[k];

    vector<int> input_I_asym_pops(n_pop, 0); //OK
    vector<int> input_I_presym_pops(n_pop, 0); //OK
    vector<int> input_I_sym_pops(n_pop, 0); //OK
    vector<int> input_I_home_pops(n_pop, 0); //OK
    vector<int> input_I_hosp_pops(n_pop, 0); //OK
    vector<int> input_I_icu1_pops(n_pop, 0); //OK
    vector<int> input_I_icu2_pops(n_pop, 0); //OK
    vector<int> input_R_pops(n_pop, 0); //OK
    vector<int> input_D_pops(n_pop, 0); //OK

    // Convert R0 to a list
    int total_windows = input_r0.size(); //not_used
    double *r0 = (double *)malloc(n_pop * total_windows * sizeof(double));
    for (int outerIndex = 0; outerIndex < n_pop; outerIndex++) {
        for (int innerIndex = 0; innerIndex < total_windows; innerIndex++) {
            r0[innerIndex + (outerIndex * total_windows)] = input_r0[innerIndex];
        }
    }

    COVID19ParamStruct params;
    // General parameters
    params.input_realz_seeds = &input_realz_seeds[0];
    params.input_census_area = &input_census_area[0];
    params.input_dist_vec = &input_dist_mat[0];
    params.input_r0 = r0;
    params.input_dist_param = &input_dist_param[0];
    params.input_m = &input_m[0];
    params.input_imm_frac = &input_imm_frac[0];
    params.input_hosp_rate = &input_hosp_rate[0];
    params.input_icu_rate = &input_icu_rate[0];
    params.input_death_rate = &input_death_rate[0];
    params.input_recov_hosp = &input_recov_hosp[0];
    params.input_window_length = &input_window_length[0];
    params.total_windows = total_windows;
    params.n_realz = n_realz;
    params.n_pop = n_pop;
    params.t_max = t_max;
    params.tau = 1;
    params.n_equations = 11;
    params.trans_type = trans_type;
    params.stoch_sd = stoch_sd;
    params.dd_trans_monod_k = dd_trans_monod_k;
    // COVID19 parameters
    params.input_N_pops = &input_N_pops[0];
    params.input_S_pops = &input_S_pops[0];
    params.input_E_pops = &input_E_pops[0];
    params.input_I_asym_pops = &input_I_asym_pops[0];
    params.input_I_presym_pops = &input_I_presym_pops[0];
    params.input_I_sym_pops = &input_I_sym_pops[0];
    params.input_I_home_pops = &input_I_home_pops[0];
    params.input_I_hosp_pops = &input_I_hosp_pops[0];
    params.input_I_icu1_pops = &input_I_icu1_pops[0];
    params.input_I_icu2_pops = &input_I_icu2_pops[0];
    params.input_R_pops = &input_R_pops[0];
    params.input_D_pops = &input_D_pops[0];
    params.delta = delta;
    params.recov_a = recov_a;
    params.recov_p = recov_p;
    params.recov_s = recov_s;
    params.recov_home = recov_home;
    params.recov_icu1 = recov_icu1;
    params.recov_icu2 = recov_icu2;
    params.asym_rate = asym_rate;
    params.sym_to_icu_rate = sym_to_icu_rate;
    params.frac_beta_asym = frac_beta_asym;
    params.frac_beta_hosp = frac_beta_hosp;

	int status = covid19_model( &params,
                                &out_pops_seed[0],
                                &out_pops_pop[0],
                                &out_pops_time[0],
                                &out_pops_S_pop[0],
                                &out_pops_E_pop[0],
                                &out_pops_I_asym_pop[0],
                                &out_pops_I_presym_pop[0],
                                &out_pops_I_sym_pop[0],
                                &out_pops_I_home_pop[0],
                                &out_pops_I_hosp_pop[0],
                                &out_pops_I_icu1_pop[0],
                                &out_pops_I_icu2_pop[0],
                                &out_pops_R_pop[0],
                                &out_pops_D_pop[0],
                                &out_events_pos[0],
                                &out_events_sym[0],
                                &out_events_total_hosp[0],
                                &out_events_total_icu[0],
                                &out_events_n_death[0],
								beta_vec);

    free(r0);

    OutputParams out;
    out.out_total_hosp = new int[nrow];
    std::copy (out_events_total_hosp.begin(), out_events_total_hosp.end(), out.out_total_hosp);
    out.out_pops_time = new int[nrow];
    std::copy (out_pops_time.begin(), out_pops_time.end(), out.out_pops_time);
    out.out_pops_seed = new int[nrow];
    std::copy (out_pops_seed.begin(), out_pops_seed.end(), out.out_pops_seed);
    out.out_pops_pop = new int[nrow];
    std::copy (out_pops_pop.begin(), out_pops_pop.end(), out.out_pops_pop);

//    printf("COVID19 status %d\n", status);
    return out;
}
