#include <stdlib.h>
#include <time.h>

void calculate_electic_field(N, phi, E);
void calculate_electron_density(N, r, ne);
void calculate_potential(N, phi, rho, kappa);
double distribution(double vb); // Sample from velocity distribution of beam
void evolve_solution(double &t, int N, double *r, double *v, void (*rhs_eval)(double &, int, double *, double *, double *, double *), double dt); // Evolve solution by one timestep
void normalize_coordinates(int N, double *r, double *v);
void rhs_eval(double &t, int N, int J, double *r, double *v, double dt); // calculate derivative for evolution proces

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

  // Now actually run it

  int num_steps = (int)floor(tmax / dt);
  for(int i = 0; i < num_steps; ++i) {
    evolve_solution(t, N, r, v, rhs_eval, dt);
    normalize_coordinates(N, r, v);
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

// Numerical integrator

void evolve_solution(double &t, int N, double *r, double *v, void (*rhs_eval)(double, int, double *, double *, double *, double *), double dt) {
  
}

void normalize_coordinates(int N, double *r, double *v) {
  for(int i = 0; i < N; ++i) {
    if(r[i] < 0.0) 
      r[i] += L;
    if(r[i] > L)
      r[i] -= L;
  }
}

// Calculate derivative

void rhs_eval(double &t, int N, int J, double *r, double *v, double *rdot, double *vdot) {
  normalize_coordinates(N, r, v);

  double ne[N];
  calculate_electron_density(N, r, ne);
  
  double n0 = (double)N/L;
  double rho[J];
  for(int j = 0; j < J; ++j)
    rho[J] = ne[j] / n0 - 1.0;
  double kappa = 2.0 * M_PI / L;
  double phi[N];
  calculate_potential(N, phi, rho, kappa);

  double E[N];
  calculate_electic_field(N, phi, E);

  for(int i = 0; i < N; ++i) {
    double dx = L / (double)J;
    int j = (int)(floor(r[i]/dx));
    double y = r[i]/dx - (double)j;
    
    double Efield = E[j] * (1.0 - y) + E[((j+1) == J) ? 0 : j + 1] * y;
    
    rdot[i] = v[i];
    vdot[i] = -Efield;
  }
  
}
