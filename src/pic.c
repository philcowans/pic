#include <fftw.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void calculate_electic_field(int J, int L, double *phi, double *E);
void calculate_electron_density(int N, int J, int L, double *r, double *ne);
void calculate_potential(int N, int J, double *phi, double *rho, double kappa);
double distribution(double vb); // Sample from velocity distribution of beam
void evolve_solution(double t, int N, int J, int L, double *r, double *v, void (*rhs_eval)(double, int, int, int, double *, double *, double *, double *), double dt); // Evolve solution by one timestep
void fft_forward(int J, double *f, double *Fr, double *Fi);
void fft_backward(int J, double *Fr, double *Fi, double *f);
void normalize_coordinates(int N, int L, double *r, double *v);
void rhs_eval(double t, int N, int J, int L, double *r, double *v, double *rdot, double *vdot); // calculate derivative for evolution proces

int main(int argc, char **argv) {
  
  // Simulation parameters

  int L = 100; // Domain length (Debye lengths)
  int N = 20000; // Number of electrons to simulate
  int J = 1000; // Number of grid points

  double vb = 3.0; // Beam velocity
  double dt = 0.1; // Time step (inverse plasma frequency units)
  double tmax = 20.0; // Time scale for simulation

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
    evolve_solution(t, N, J, L, r, v, rhs_eval, dt);
    normalize_coordinates(N, L, r, v);
  }

  for(int i = 0; i < N; ++i) {
    printf("%f, %f\n", r[i], v[i]);
  }

  return 0;
}

void calculate_electic_field(int J, int L, double *phi, double *E) {
  double dx = L / (double)J;

  for(int j = 1; j < J-1; ++j)
    E[j] = ((phi[j-1] - phi[j+1]) / 2.0) / dx;

  E[0] = ((phi[J-1] - phi[1]) / 2.0) / dx;
  E[J-1] = ((phi[J-2] - phi[0]) / 2.0) / dx;
}

void calculate_electron_density(int N, int J, int L, double *r, double *ne) {
  double dx = L / (double)J;
  
  for(int j = 0; j < J; ++j) {
    ne[j] = 0.0;
  }

  for(int i = 0; i < N; ++i) {
    int j = (int)floor(r[i] / dx);
    double y = r[i] / dx - (double)j;
    ne[j] += (1.0 - y) / dx;
    if(j+1 == J) 
      ne[0] += y / dx;
    else 
      ne[j+1] += y / dx;
  }
}

void calculate_potential(int N, int J, double *phi, double *rho, double kappa) {
  double Vr[J], Vi[J], Ur[J], Ui[J];

  // Fourier transform source term
  fft_forward(J, rho, Vr, Vi);

  Ur[0] = Ui[0] = 0.0;

  for(int j = 1; j <= J/2; ++j) {
    Ur[j] = -Vr[j] / (double)(j * j) / kappa / kappa;
    Ui[j] = -Vi[j] / (double)(j * j) / kappa / kappa;
  } 

  for (int j = J/2; j < J; ++j) {
    Ur[j] = Ur[J-j];
    Ui[j] = -Ui[J-j];
  }

  fft_backward(J, Ur, Ui, phi);
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

void evolve_solution(double t, int N, int J, int L, double *r, double *v, void (*rhs_eval)(double, int, int, int, double *, double *, double *, double *), double dt) {
  double rk1[N];
  double rk2[N];
  double rk3[N];
  double rk4[N];
  double rf[N];
  double rdot[N];

  double vk1[N];
  double vk2[N];
  double vk3[N];
  double vk4[N];
  double vf[N];
  double vdot[N];

  // Zeroth intermediate step 
  (*rhs_eval)(t, N, J, L, r, v, rdot, vdot);
  for(int j = 0; j < N; ++j) {
    rk1[j] = dt * rdot[j];
    vk1[j] = dt * vdot[j];
    rf[J] = r[J] + rk1[J] / 2.0;
    vf[J] = v[J] + vk1[J] / 2.0;
  }

  // First intermediate step 
  (*rhs_eval)(t + dt/2.0, N, J, L, rf, vf, rdot, vdot);
  for(int j = 0; j < N; ++j) {
    rk2[j] = dt * rdot[j];
    vk2[j] = dt * vdot[j];
    rf[J] = r[J] + rk2[J] / 2.0;
    vf[J] = v[J] + vk2[J] / 2.0;
  }

  // Second intermediate step
  (*rhs_eval)(t + dt/2.0, N, J, L, rf, vf, rdot, vdot);
  for(int j = 0; j < N; ++j) {
    rk3[j] = dt * rdot[j];
    vk3[j] = dt * vdot[j];
    rf[J] = r[J] + rk3[J];
    vf[J] = v[J] + vk3[J];
  }

  // Third intermediate step 
  (*rhs_eval)(t + dt, N, J, L, rf, vf, rdot, vdot);
  for(int j = 0; j < N; ++j) {
    rk4[j] = dt * rdot[j];
    vk4[j] = dt * vdot[j];
  }

  // Actual step 
  for(int j = 0; j < N; ++j) {
    r[j] += rk1[j] / 6.0 + rk2[j] / 3.0 + rk3[j] / 3.0 + rk4[j] / 6.0;
    v[j] += vk1[j] / 6.0 + vk2[j] / 3.0 + vk3[j] / 3.0 + vk4[j] / 6.0;
  }

  // TODO: This isn't actually returned (not an issue here)
  t += dt;
}

// Calculates Fourier transform of array f in arrays Fr and Fi
void fft_forward(int J, double *f, double *Fr, double *Fi) {
  fftw_complex ff[J], FF[J];

  for (int j = 0; j < J; ++j) {
    ff[j].re = f[j]; 
    ff[j].im = 0.0;
  }

  fftw_plan p = fftw_create_plan(J, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_one(p, ff, FF);
  fftw_destroy_plan(p); 

  for (int j = 0; j < J; ++j) {
    Fr[j] = FF[j].re; 
    Fi[j] = FF[j].im;
  }

  // Normalize data
  for (int j = 0; j < J; ++j) {
    Fr[j] /= (double)J;
    Fi[j] /= (double)J;
  }
}

// Calculates inverse Fourier transform of arrays Fr and Fi in array f
void fft_backward(int J, double *Fr, double *Fi, double *f) {
  fftw_complex ff[J], FF[J];

  for(int j = 0; j < J; ++j) {
    FF[j].re = Fr[j]; 
    FF[j].im = Fi[j];
  }

  fftw_plan p = fftw_create_plan(J, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_one(p, FF, ff);
  fftw_destroy_plan(p); 

  for(int j = 0; j < J; ++j)
    f[j] = ff[j].re; 
}

void normalize_coordinates(int N, int L, double *r, double *v) {
  for(int i = 0; i < N; ++i) {
    if(r[i] < 0.0) 
      r[i] += L;
    if(r[i] > L)
      r[i] -= L;
  }
}

// Calculate derivative

void rhs_eval(double t, int N, int J, int L, double *r, double *v, double *rdot, double *vdot) {
  normalize_coordinates(N, L, r, v);

  double ne[N];
  calculate_electron_density(N, J, L, r, ne);
  
  double n0 = N/(double)L;
  double rho[J];
  for(int j = 0; j < J; ++j)
    rho[j] = ne[j] / n0 - 1.0;
  double kappa = 2.0 * M_PI / L;
  double phi[J];
  calculate_potential(N, J, phi, rho, kappa);

  double E[J];
  calculate_electic_field(J, L, phi, E);

  for(int i = 0; i < N; ++i) {
    double dx = L / (double)J;
    int j = (int)(floor(r[i]/dx));
    double y = r[i]/dx - (double)j;
    
    double Efield = E[j] * (1.0 - y) + E[((j+1) == J) ? 0 : j + 1] * y;
    
    rdot[i] = v[i];
    vdot[i] = -Efield;
  }
}
