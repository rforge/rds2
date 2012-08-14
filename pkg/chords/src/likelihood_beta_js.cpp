#include <vector>
#include <limits>
#include <iostream>
#include <math.h>
#include <R.h>

using namespace std;

extern "C"{
	void likelihood_beta_js(  int *sample,
					  int *Sij,
					  int *S,
					  double *beta_js,
					  double *Nj,
					  int *observed_degrees,
					  int *n,
					  int *N,
					  int *N_observed,
					  bool *arc,
					  double *result ) {				
		
		bool			first = true;
		
		*result=0.0;
		int temp_degree=0;

		 for(int i=0; i < *n; ++i)  { // go over observations
			for(int j=0; j < *N_observed; ++j) { // go over observed degrees
				temp_degree = observed_degrees[j];
								*result += (sample[i] == temp_degree ) ?
									(log(beta_js[j]) + log((double)S[i]) + log(Nj[temp_degree-1] - Sij[j + i*(*N_observed) ]) ) :
									(log(1.0 - beta_js[j] * S[i] * ( Nj[temp_degree-1] - Sij[j + i*( *N_observed)] )));
			}
			if(sample[i]>0 && first) 
				first = false;
		 }
		
		if(isnan(*result))
			*result = -numeric_limits<double>::infinity();			
	}
}
