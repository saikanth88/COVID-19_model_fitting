# COVID-19 parameter estimation

## About the Project
In this project, we propose a method to estimate time-varying ransmission rates and host movement rates in a stochastic and spatially-explicit model of SARS-CoV-2 using  grid search algorithm with sliding window technique. You can find the SARS-CoV-2 model used for the implementation of the code in [medrxiv](https://www.medrxiv.org/content/10.1101/2021.05.13.21256216v1). 

## Getting Started
All the codes were written using C++ and R languages. The main codes of parameter estimation and computing new hospitalizations were written using C++. The R code is used to generate figures of estimated parameters and hospitalizations.

## Usage

## Roadmap
### ./data
* [region]_R0_hosp_count.csv - daily new hospitalizations
* [region]_countyData.csv - spatial data (population size, population density, area, latitude, longitude, etc.)
* [region]_daily_beta.csv, region_beta_vals.csv - daily transmission rate (beta) values
* region = `low` for rural, region = `mid` for suburban, and region = `high` for urban

### ./parameter_estimation
* Select the case: estimate only transmission rate `(beta)`, estimate only movement rate `(m)`, or estimate both transmission rate and movement rate `(beta_m)`.
* grid_search_covid_[case].cpp contains the main C++ file.
* Set the `METHOD` to 0 or 1 to choose the likelihood calculation method.
* Set number of realizations `(n_realz)` and number of grid searches `(n_searches)`.
* Set `const string region` to `low`, `mid`, or `high` to select the region.
* Include [region]_R0_hosp_count.csv, [region]_countyData.csv, and [region]_daily_beta.csv (only need to estimate m) files in the same folder.
* Use the `Makefile` to compile and run the code.
* The code will generate "final_output.csv" with weekly parameter estimations.

Note: *You may need to uncomment the line `#include <uuid/uuid.h>` in county_param.h and the lines related to creating a filename using `uuid` in grid_search_covid_[case].cpp to run the code in parallel.*
  
### ./model_hospitalizations
* Concatenate all outputs into one .csv file (e.g.: [region]_full_covid_list_5000_all_pops.csv), if the parameter estimation was implemented in parallel.
* Select the case: model hospitalizations using estimated beta, hospitalizations using estimated m, or hospitalizations using estimated beta and m.
* covid_hospitalizations_[case].cpp contains the main C++ file.
* Include [region]_countyData.csv, and [region]_daily_beta.csv (only need to estimate m) files in the same folder.
* Change the filename covid_hospitalizations_[case] in `Makefile` according to the case.
* Compile and run the code.
* The code will generate "params_all_counts_5000_all_pops.csv" with daily new hospitalizations, total hospitalizations,and daily new deaths.

### ./output_generation
* User must have output files for all three regions for a specific case to run this code.
* Select the R file with the correct case: generate_output_[case].R
* Set the path to the .csv files in the code.
* Files need to generate the figures: [region]_R0_hosp_count.csv, [region]_countyData.csv, [region]_beta_vals.csv, [region]_full_covid_list_5000_all_pops.csv, [region]_params_all_counts_5000_all_pops.csv
* The code assumes both [region]_full_covid_list_5000_all_pops.csv and [region]_params_all_counts_5000_all_pops.csv files are stored in a sub-folder in the main folder with region specific .csv files.
* The code will create plots of weekly transmission rate estimations, weekly movement rate esitmations, and new hospitalizations for the three regions. 

## License

## Contact
Saikanth Ratnavale: saikanth.ratnavale@nau.edu, saik8801@gmail.com

Joseph Mihaljevic: joseph.mihaljevic@nau.edu

## Acknowledgments
