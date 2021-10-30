debug:
	g++ -I/packages/gsl/2.6/include -O3 -c -g covid_hospitalizations_beta.cpp
	g++ -L/packages/gsl/2.6/lib covid_hospitalizations_beta.o -luuid -lgsl -lgslcblas -lm -o covid_hospitalizations_beta

clean:
	rm -f *.o
	rm -f *.exe
	rm -f FixMCInt
