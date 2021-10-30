//Calculating county specific parameters for covid19_driver()

using namespace std;

double distance(double lat1, double long1, double lat2, double long2);
double *dist_matrix(double *ilon, double *ilat, int count);

CountyParams county_data(const string region){

	CountyParams params;

	int row = 0;

	ifstream fin;

	string county_file = region+"_countyData.csv";
	string line;
	fin.open(county_file);
	fin.ignore(1000, '\n');
	//ignore the first 1000 characters, or until first \n, whichever is met first
	while(getline(fin,line))
		row++;

	fin.close();

	printf("rows: %d\n", row);

	int n_pop = row;

	string	pop_size, pop_density, pop_ID, area, type, lon, lat, county;
	double* ipop_size = new double[row];
	double* ipop_density = new double[row];
	double* iarea = new double[row];
	double* ilon = new double[row];
	double* ilat = new double[row];

	int* round_pop = new int[row];


	fin.open(county_file);
	fin.ignore(1000, '\n');

	row = 0;
	while(row < n_pop)
	{
		getline( fin, pop_size, ',' );
		getline( fin, pop_density, ',' );
		getline( fin, pop_ID, ',' );
		getline( fin, area, ',' );
		getline( fin, type, ',' );
		getline( fin, lon, ',' );
		getline( fin, lat, ',' );
		getline( fin, county, '\n' );
		stringstream( pop_size ) >> ipop_size[row];
		stringstream( pop_density ) >> ipop_density[row];
		stringstream( area ) >> iarea[row];
		stringstream( lon ) >> ilon[row];
		stringstream( lat ) >> ilat[row];

		round_pop[row] = (int)(round(ipop_size[row]));
//      printf("%.3f %.3f %.3f %.3f %.3f\n", ipop_size[row], ipop_density[row], iarea[row], ilon[row], ilat[row]);
//		printf("%f %d\n", ipop_size[row],round_pop[row]);
		row++;
	}

	fin.close();

//	double **distm = dist_matrix(ilon, ilat, n_pop);
	double *distv = dist_matrix(ilon, ilat, n_pop);

///////////////////////////////////////////////////////////////////////////
	int max_pop = *std::max_element(round_pop, round_pop + n_pop);
	int imax;
	double mlat, mlon;
	for (int i = 0; i < n_pop; i++){
        if (round_pop[i] == max_pop)
            imax = i;
	}
	mlat = ilat[imax], mlon = ilon[imax];

//    cout<<"max_pop:"<<max_pop<<" mlat:"<<mlat<<" mlon:"<<mlon<<"\n";

	double *dist_diff = new double[n_pop];
	int *dist_ID = new int[n_pop];
    for (int i = 0; i < n_pop; i++){
        dist_diff[i] = distance(ilat[i], ilon[i], mlat, mlon);
        if (dist_diff[i] == 0)
            dist_ID[i] = 1;
        else if ((dist_diff[i] > 0)&&(dist_diff[i] <= 50))
            dist_ID[i] = 2;
        else if ((dist_diff[i] > 50)&&(dist_diff[i] <= 100))
            dist_ID[i] = 3;
        else if ((dist_diff[i] > 100)&&(dist_diff[i] <= 150))
            dist_ID[i] = 4;
        else if ((dist_diff[i] > 150)&&(dist_diff[i] <= 200))
            dist_ID[i] = 5;
        else if ((dist_diff[i] > 200)&&(dist_diff[i] <= 250))
            dist_ID[i] = 6;
        else
            dist_ID[i] = 7;
//        cout<<dist_diff[i]<<"\t"<<dist_ID[i]<<"\n";
	}

	std::vector<int> E_pops_vec(n_pop, 0);

    if(region == "low")
        E_pops_vec[imax] = 15;
    else if(region == "mid")
        E_pops_vec[imax] = 20;
    else if(region == "high")
        E_pops_vec[imax] = 25;

///////////////////////////////////////////////////////////////////////////

	params.n_pop 		= n_pop;
	params.pop_N  		= new int[n_pop];
	params.E_pops 		= new int[n_pop];
	params.census_area	= new double[n_pop];
	params.dist_vec		= new double[n_pop*n_pop];
	params.dist_ID      = new int[n_pop];

	memcpy(params.pop_N, round_pop, n_pop*sizeof(int));
//	memcpy(params.E_pops, E_pops, n_pop*sizeof(int));
	std::copy(E_pops_vec.begin(), E_pops_vec.end(), params.E_pops);
	memcpy(params.census_area, iarea, n_pop*sizeof(double));
	memcpy(params.dist_vec, distv, n_pop*n_pop*sizeof(double));
	memcpy(params.dist_ID, dist_ID, n_pop*sizeof(int));

//	for (int k = 0; k < n_pop; k++){
//		printf("E_pops[%d]: %d\n", k,params.E_pops[k]);
//	}

	delete[] ipop_size;
	delete[] ipop_density;
	delete[] iarea;
	delete[] ilon;
	delete[] ilat;
	delete[] distv;
	delete[] round_pop;
    delete[] dist_diff;
	delete[] dist_ID;

	return params;

	delete[] params.pop_N;
	delete[] params.E_pops;
	delete[] params.census_area;
	delete[] params.dist_vec;
    delete[] params.dist_ID;
}

//https://www.geeksforgeeks.org/program-distance-two-points-earth/
// converting degrees to radians
double toRadians(const double degree)
{
	double one_deg = (M_PI)/(180.0);
	return (one_deg * degree);
}

double distance(double lat1, double long1, double lat2, double long2)
{
	// Convert the latitudes and longitudes from degree to radians.
	lat1  = toRadians(lat1);
	long1 = toRadians(long1);
	lat2  = toRadians(lat2);
	long2 = toRadians(long2);

	// Haversine Formula
	double dlong = long2 - long1;
	double dlat  = lat2 - lat1;

	double ans = pow(sin(dlat/2),2) + cos(lat1) * cos(lat2) * pow(sin(dlong/2),2);

	ans = 2 * asin(sqrt(ans));

	// Radius of Earth in Kilometers, R = 6371, Use R = 3956 for miles
	double R = 6378.137;

	// Calculate the result
	ans = ans * R;

	return ans;
}

double *dist_matrix(double* ilon, double* ilat, int n_pop)
{
	int rows = n_pop;
	int cols = n_pop;
	double **dist_mat = new double *[rows];
	double *dist_vec = new double[rows*cols];

	for(int i = 0; i < rows; i++)
		dist_mat[i] =  new double [cols];

//https://cran.r-project.org/web/packages/geosphere/geosphere.pdf
//https://rdrr.io/cran/geosphere/man/distm.html
//use distGeo, current is Distcosine

	for (int i = 0; i<rows; i++){
		for (int j = 0; j<cols; j++){
			dist_mat[i][j] = distance(ilat[i], ilon[i], ilat[j], ilon[j]);
			dist_vec[i + j*n_pop] = dist_mat[i][j];
		}
	}

	return dist_vec;

	for (int i = 0; i < rows; i++)
		delete[] dist_mat[i];

	delete[] dist_mat;
	delete[] dist_vec;
}
