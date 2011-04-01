//     compile with: R CMD SHLIB rw.c

#include <R.h>
#include <Rmath.h>

void work(int *Bin,
	  int *Pin,
	  int *nin, 
	  double *f1in,
	  int *m, 
	  int *m_partial_sum,
	  int *m_sumin,
	  int *m_maxin,
	  int *M,
	  double *M_barin, 
	  double *lam1in,  
	  double *y_bar,
	  double *y,
	  double *tilde_bar)
{
    GetRNGstate();

    // declare stuff that already exists
    int n = nin[0]; // number of sampled PSUs
    double M_bar = M_barin[0]; // average number of SSUs per PSU
    double lam1 = lam1in[0]; // sqrt(n/(n-1) * (1-n/N))
    double f1 = f1in[0]; // n/N
    int P = Pin[0];  // number of variables
    int m_sum = m_sumin[0]; // number of SSUs in whole sample
    int m_max = m_maxin[0]; // largest individual secondary sample
    int B = Bin[0]; // how many bootstraps to do

    // declare stuff that doesn't exist yet
    int star_psu; // resampled PSU indices
    int star_m; // number of sampled SSUs in resampled PSUs
    int star_M; // number of total SSUs in resampled PSUs
    double star_y[m_max]; // only first star_m will be filled
    double tilde_y[m_max]; // only first star_m will be filled 
    double f2; // secondary sampling fraction
    double lam2; // function of sampling fractions and sizes
    int v[m_max]; // randomly selected indices of SSUs
    double star_y_sum; // sum of a few star_y's
    double syt; // estimated (HT) total of a few star_y's

    // bootstrap loop will be below this line
    for(int b = 0; b < B; b++){
      for(int i = 0; i < n; i++){
	// figure out which PSU we're working with and set design variables
	double u = unif_rand();
	star_psu = floor(u*n);
	star_M = M[star_psu];
	star_m = m[star_psu];
	f2 = ((double) star_m)/((double) star_M);
	lam2 = sqrt(f1 * (1-f2) * ((double) star_m)/((double) star_m-1)); // weird fct	  

	// select SSU indices
	for(int j = 0 ; j < star_m; j++){
	  u = unif_rand();
	  v[j] = floor(u*((double) star_m));
	}
	
	for(int p = 0; p < P; p++){
	  // for the pth variable, find star_y and star_y_total
	  star_y_sum = 0; // sum of star_y
	  for(int j = 0; j < star_m; j++){
	    star_y[j] =  y[p*m_sum + m_partial_sum[star_psu]+v[j]];
	    star_y_sum = star_y_sum + star_y[j];
	  }
	  syt = star_y_sum*((double) star_M)/((double) star_m); // (HT) estimated PSU total
	  
	  // for the pth variable, find tilde_y and add up contribution to tilde_bar
	  for(int j = 0 ; j < star_m; j++){
	    tilde_y[j] = 
	      y_bar[p] + lam1*(syt/M_bar - y_bar[p]) + 
	      lam2*(((double) star_M)*star_y[j]/M_bar - syt/M_bar);
	    tilde_bar[b*P+p] = tilde_bar[b*P+p] + (tilde_y[j]/((double) star_m))/((double) n);
	  }
	}
      }
      PutRNGstate();
    }
}
