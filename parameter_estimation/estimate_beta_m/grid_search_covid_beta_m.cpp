#include "county_param.h"
#include "county_data.cpp"
#include "covid19_driver.cpp"

using namespace std;

//************* MAKE THE FOLLOWING CHANGES BEFORE YOU RUN *****************

//************** Sets the likelihood calculation method *******************
// 1 - with pop groups, 0 - without pop groups (default)
#define METHOD 1

//******************** Sets the region of interest ************************
// low, mid, or high
const string region = "low";

//*************************************************************************

gsl_rng *rand2;

CountyParams par = county_data(region);
int n_pop  = par.n_pop;

double LHood_Poisson(int **data_hosp, int n_days, int win_count, int window_size, double *param, BestParams best)
{
  int t_max = 7*win_count + 7*window_size; //in days(previous + current)
  int start_day, end_day;

  if((n_days - t_max) < 7)
    t_max = t_max + (n_days - t_max); //add any remaining days to the last window

  double *beta_vec = new double[t_max];
  double *mgrt_vec = new double[t_max];

  if(param != NULL){ //to calculate best params
    if (win_count > 0){
      for(int i = 0; i < win_count; i++){ //previous windows
        for(int j = 7*i; j < 7*(i + 1); j++){
          beta_vec[j] = best.beta_best[i];
          mgrt_vec[j] = best.mgrt_best[i];
        }
      }
    }

    for(int j = 7*win_count; j < t_max; j++){ //current 3 week window
      beta_vec[j] = param[0];
      mgrt_vec[j] = param[1];
    }
    start_day = 7*win_count;
    end_day   = t_max;
  }else{ //to calculate full_LHood
    for(int i = 0; i <= win_count; i++){
      for(int j = 7*i; j < 7*(i + 1); j++){
        beta_vec[j] = best.beta_best[i];
        mgrt_vec[j] = best.mgrt_best[i];
      }
    }
    start_day = 0;
    end_day   = 7*(win_count + 1); //all windows computed before
  }

  int n_realz = 2;

  OutputParams model = covid19_driver(end_day, n_realz, beta_vec, mgrt_vec, par);

  int real_counter = 0;
  int day_counter  = 0;
  int pop_counter  = 0;

  int county_hosp[n_realz][end_day][n_pop];

  for(int i = 0; i < n_realz*n_pop*end_day; i++){
    real_counter = model.out_pops_seed[i] - 1; //(-1) change the starting index from 1 to 0
    day_counter  = model.out_pops_time[i] - 1;
    pop_counter  = model.out_pops_pop[i] - 1;

    county_hosp[real_counter][day_counter][pop_counter] = model.out_total_hosp[i];
  }

  int k_max = n_pop;

  //for METHOD 1
  int id_max = *std::max_element(par.dist_ID, par.dist_ID + n_pop);
  int data_hosp2[n_days][id_max];
  int county_hosp2[n_realz][end_day][id_max];

  if (METHOD == 1){

    k_max = id_max;

    std::fill_n(&data_hosp2[0][0], n_days*id_max, 0); //initialize all elements to zero
    std::fill_n(&county_hosp2[0][0][0], n_realz*end_day*id_max, 0);

    for (int i = 0; i < n_days; i++){
      for(int j = 0; j < n_pop; j++){
        data_hosp2[i][par.dist_ID[j] - 1] += data_hosp[i][j]; //shift the group to the left by 1
      }
    }


    for (int i = 0; i < n_realz; i++) {
      for(int j = 0; j < end_day; j++){
        for(int k = 0; k < n_pop; k++){
          county_hosp2[i][j][par.dist_ID[k] - 1] += county_hosp[i][j][k];
        }
      }
    }

  } //end method 1

  double lhood;

  double *sum_log_lhood = new double[n_realz];
  double avg_log_lhood = 0.0;

  memset(sum_log_lhood, 0.0, n_realz * sizeof(double));

  for (int i = 0; i < n_realz; i++) { //check likelihood only for the current window
    for(int j = start_day; j < end_day; j++){
      for (int k = 0; k < k_max; k++) {

        if(METHOD == 1)
            lhood = gsl_ran_poisson_pdf(data_hosp2[j][k], county_hosp2[i][j][k]+0.001);
        else
            lhood = gsl_ran_poisson_pdf(data_hosp[j][k], county_hosp[i][j][k]+0.001);

        if((isinf(log(lhood))!=0) || (isnan(log(lhood))!=0)){
          // If the lhood is so bad that it gives NaN or -INF:
          sum_log_lhood[i] += -700.0;
        }else{
          sum_log_lhood[i] += log(lhood);
        }

      } //end pop or group(k) loop
    } //end day(j) loop
  } //end realize(i) loop

  for (int j = 0; j < n_realz; j++) {
    avg_log_lhood = avg_log_lhood + sum_log_lhood[j];
  }

  avg_log_lhood = avg_log_lhood/n_realz;

  delete[] beta_vec;
  delete[] mgrt_vec;
  delete[] sum_log_lhood;

  return(avg_log_lhood);
}


int main()
{
  BestParams best;

  int n_rows = 0;

  //////////////////////// READ CSV FILE ////////////////////////
  ifstream fin;
  //Extract rows of one county
  string data_file = region+"_R0_hosp_count.csv";
  string line;
  fin.open(data_file);
  fin.ignore(1000, '\n');
  while(getline(fin,line))
    n_rows++;
  fin.close();

  int n_days = n_rows/n_pop;
  printf("n_days: %d\n", n_days);

  string str_day, str_pop, str_count;
  int day, pop;

  int **n_hosp = new int *[n_days];

  for(int i = 0; i < n_days; i++)
	n_hosp[i] =  new int[n_pop];

  fin.open(data_file);
  fin.ignore(1000, '\n');
  for(int i = 0; i < n_rows; i++){
    getline( fin, line);
    stringstream iss(line);
    getline(iss, str_day, ',');
    getline(iss, str_pop, ',');
    getline(iss, str_count, '\n');
    stringstream(str_day) >> day;
    stringstream(str_pop) >> pop;
    stringstream(str_count) >> n_hosp[day - 1][pop - 1];
//    cout<<day<<"\t"<<pop<<"\t"<<n_hosp[day - 1][pop - 1]<<"\n";
  }
  fin.close();

  ///////////////////////////////////////////////////////////////

  const gsl_rng_type* T1 = gsl_rng_default;
  rand2 = gsl_rng_alloc(T1);
  // Set the seed:
  // (based on time and process id)
  gsl_rng_set(rand2, (time(0)+getpid()));

  int n_param = 2;

  double *param_min  = new double[n_param];
  double *param_max  = new double[n_param];
  double *step_size  = new double[n_param];
  double *param_best = new double[n_param];
  double *param_temp = new double[n_param];

  param_min[0] = 0.001, param_max[0] = 0.6; //beta
  param_min[1] = 0.0,   param_max[1] = 0.08; //m

  int window_size = 3; //in weeks
  int n_weeks = (int)(n_days/7.0);
  int n_windows = n_weeks - window_size + 1;
  best.beta_best = new double[n_windows];
  best.mgrt_best = new double[n_windows];
  memset(best.beta_best, 0.0, n_windows*sizeof(double));
  memset(best.mgrt_best, 0.0, n_windows*sizeof(double));

//  ofstream fout1;
//  fout1.open ("test_param_output.csv");

  for(int win_count = 0; win_count < n_windows; win_count++){
    if (win_count == 0){
      for(int k = 0; k < n_param; k++){
        // This holds temporary parameter values
        param_temp[k] = param_min[k] + (param_max[k] - param_min[k])/ gsl_ran_flat(rand2, 1, 10);
      }
    }else{
        param_temp[0] = best.beta_best[win_count-1];
        param_temp[1] = best.mgrt_best[win_count-1];
    }

    double adder_1, adder_2, n_steps;
    double temp_param;
    double rand_param;
    double temp_max = 0.0;
    double max_lhood = -1000000.0;

    int PC, starter;
    int n_searches = 1;

    temp_max = LHood_Poisson(n_hosp, n_days, win_count, window_size, param_temp, best);
    max_lhood = temp_max;
    // Copy all the temp params to the param_best, so it isn't blank
    for(int k = 0; k < n_param; k++){
      memcpy(&param_best[k], &param_temp[k], sizeof(double));
    }

    for(int search_ID = 0; search_ID < n_searches; search_ID++){
      adder_1 = gsl_rng_uniform(rand2);
      adder_2 = gsl_rng_uniform(rand2);
      n_steps = 5*adder_1 + 5*adder_2;

      for(int k = 0; k < n_param; k++){
        step_size[k] = (param_max[k] - param_min[k])/n_steps;
      }

      starter = gsl_rng_uniform_int(rand2, n_param);

      // Cycle through the parameters
      for(int index = starter;  index < (n_param + starter);  index++){
        PC = index % n_param;

        for(temp_param = param_min[PC]; temp_param < param_max[PC]; temp_param += step_size[PC]){
          if((temp_param + step_size[PC]) < param_max[PC])
             rand_param = gsl_ran_flat(rand2, temp_param, temp_param + step_size[PC]);
          else
             rand_param = gsl_ran_flat(rand2, temp_param, param_max[PC]);

          param_temp[PC] = rand_param;
          temp_max = LHood_Poisson(n_hosp, n_days, win_count, window_size, param_temp, best);

          if(temp_max >= max_lhood){
            // If the new parameter is better:
            //    	    				 printf("Found a new best param value.\n");
            //    	    				 printf("New param_temp[%d]: %f\n",PC,param_temp[PC]);
            //    	    				 New max_lhood:
            max_lhood = temp_max;

            // Copy only new paramter value to param_best
            memcpy(&param_best[PC], &param_temp[PC], sizeof(double));
          }

          // Check values of arrays:
          // for(int k = 0; k < n_param; k++){
          //     printf("param_best[%d]: %f\n",k,param_best[k]);
          // }
        }

        //reset param_temp to current param_best
        for(int k = 0; k < n_param; k++){
          memcpy(&param_temp[k], &param_best[k], sizeof(double));
        }
      }
//      for(int k = 0; k < n_param; k++)
//        fout1 <<param_best[k]<<",";
//      fout1 << max_lhood <<"\n";
    }
    best.beta_best[win_count] = param_best[0];
    best.mgrt_best[win_count] = param_best[1];
//    cout<<beta_best[win_count]<<"\n";
    printf("win_count: %d\n", win_count);
  } //end (win_count) loop

//  fout1.close();

  //calculate the full likelihood for all n_windows considered
  double full_lhood = LHood_Poisson(n_hosp, n_days, n_windows-1, window_size, NULL, best);

//      uuid_t uuid;
//      uuid_generate_random(uuid);
//      char uuid_str[37];
//      uuid_unparse_lower(uuid, uuid_str);
//      //printf("uuid: %s\n", uuid_str);
//      string uuid_s = uuid_str;
//      string filename = uuid_s+".csv";

  //save the best parameter set and max log-likelihood
  ofstream fout2;
  fout2.open ("final_output.csv");//(filename);
  for(int k = 0; k < n_windows; k++)
    fout2 << best.beta_best[k] <<",";
  for(int k = 0; k < n_windows; k++)
    fout2 << best.mgrt_best[k] <<",";
  fout2 << full_lhood <<"\n";
  fout2.close();

  for (int i = 0; i < n_days; i++)
	delete[] n_hosp[i];

  delete[] n_hosp;
  delete[] param_min;
  delete[] param_max;
  delete[] step_size;
  delete[] param_best;
  delete[] param_temp;
  delete[] best.beta_best;
  delete[] best.mgrt_best;
  return 0;
}
