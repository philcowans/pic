#include <stdlib.h>
#include <time.h>

double distribution(double vb); // Sample from velocity distribution of beam

int main(int argc, char **argv) {
  
  // Simulation parameters

  int L = 100; // Domain length (Debye lengths)
  int N = 20000; // Number of electrons to simulate
  int J = 1000; // Number of grid points

  double vb = 3.0; // Beam velocity
  double dt = 0.1; // Time step (inverse plasma frequency units)
  double tmax = 10; // Time scale for simulation

  // Initialize simulation

  double t = 0.0;
  int seed = time(NULL);
  srand(seed);
  
  double r[N];
  double v[N];
  for(int i = 0; i < N; ++i) {
    r[i] = L * ((double)rand()) / ((double)RAND_MAX);
    v[i] = distribution(vb);
  }

  return 0;
}

// Rejection sampler for Maxwellian beam distributions

double distribution(double vb) {
  double fmax = 0.5 * (1.0 + exp(-2.0 * vb * vb));
  double vmin = -5.0 * vb;
  double vmax = 5.0 * vb;
  double v = vmin + (vmax - vmin) * ((double)rand()) / ((double)RAND_MAX);

  // Accept/reject value
  double f = 0.5 * (exp(-(v - vb) * (v - vb) / 2.0) +
		    exp(-(v + vb) * (v + vb) / 2.0));
  double x = fmax * ((double)rand()) / ((double)RAND_MAX);

  if (x > f) 
    return distribution (vb);
  else 
    return v;
}
