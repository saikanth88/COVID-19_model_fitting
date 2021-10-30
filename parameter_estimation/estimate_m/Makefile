debug:
	g++ -I/packages/gsl/2.6/include -O3 -c -g grid_search_covid_m.cpp
	g++ -L/packages/gsl/2.6/lib grid_search_covid_m.o -luuid -lgsl -lgslcblas -lm -o grid_search_covid_m

clean:
	rm -f *.o
	rm -f *.exe
	rm -f FixMCInt