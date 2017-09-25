#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

#define HALO_MAX 1.0E+16
#define PI 3.14159

float g7_ngal,
  g7_msub;

/*external functions
 */
float qromo(float (*func)(float), float a, float b,
	     float (*choose)(float(*)(float), float, float, int));
float midpnt(float (*func)(float), float a, float b, int n);
void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);
void sort2(int n, float arr[], int id[]);
float zbrent(float (*func)(float), float x1, float x2, float tol);
float gammln(float xx);

/* SHAM functions
 */
float func_ngal(float mag);
float func_nhalo(float m);
float subhalo_abundance(float m);
float subhalo_abundance2(float m);
float halo_abundance(float m);
float halo_abundance2(float m);
float func_subhalo1(float mhost);
float func_match_nhalo(float mmin);
float func_match_ngal(float mass);

/* local functions
 */
double poisson_prob(int n, double nave);
double continuous_poisson_prob(double n, double nbar);
float mass2magnitude(float mass, float *mag_all, float *mass_all, int *indx, int ngal);
double func_nsat(double mass, double m1, double mcut, double alpha);

float most_probable_mass(float mag, float m0, float nsat, float *mag_all, float *mass_all, int *indx, int ngal)
{
  int i,j,k,nmass,n;
  double mlo, mhi, dlogm, sigma_mag, prob_lum, prob_n, p, pmax, mmax, mass, nhalo, mag_bar;

  double m1, alpha=1, mcut, nbar;

  sigma_mag = 0.15*2.5;

  m1 = pow(10.0,12.85);
  mcut = pow(10.0,11.6);

  //return pow(nsat,1.0/alpha)*m1;

  mlo = mcut/10;
  mhi = m0*30;
  nmass = 100;
  dlogm = log(mhi/mlo)/nmass;

  pmax = 0;
  for(i=1;i<=nmass;++i)
    {
      mass = exp((i+0.5)*dlogm )*mlo;
      nhalo = halo_abundance(mass)*mass;
      mag_bar = mass2magnitude(mass,mag_all,mass_all,indx,ngal);
      prob_lum = exp(-(mag_bar-mag)*(mag_bar-mag)/(2*sigma_mag*sigma_mag));
      prob_n = continuous_poisson_prob(nsat,func_nsat(mass,m1,mcut,alpha));
      if(isnan(prob_n))
	{
	  nbar = func_nsat(mass,m1,mcut,alpha);
	  prob_n = 1.0/sqrt(2.*PI*nbar)*exp(-(nsat-nbar)*(nsat-nbar)/(2*nbar));
	}
      p = nhalo*prob_lum*prob_n;
      //p = prob_n;
      if(p>pmax) { pmax = p; mmax = mass; }
      printf("%e %e %e %e %e %e %e %e %e %e\n",m0,mag,nsat,mass,mag_bar,func_nsat(mass,m1,mcut,alpha),p,nhalo,prob_lum,prob_n);
    }
  printf("%e %e %e %e %e %e\n",m0,mmax,mmax/m0,nsat,func_nsat(mmax,m1,mcut,alpha),func_nsat(m0,m1,mcut,alpha));
  fflush(stdout);
  return mmax;
}

double func_nsat(double mass, double m1, double mcut, double alpha)
{
  if(mass<mcut)return 0;
  return pow((mass-mcut)/m1,alpha);
}

double continuous_poisson_prob(double n, double nbar)
{
  return pow(nbar,n)*exp(-nbar)/exp(gammln(n+1));
}

float mass2magnitude(float mass, float *mag_all, float *mass_all, int *indx, int ngal)
{
  int i, i1;
  
  for(i1=1;i1<=ngal;++i1)
    {
      i = indx[i1];
      if(mass>mass_all[i])return mag_all[i1];
    }
  return 0;      
}

void scatter_test()
{
  int i,j,k,n,nmax,nmass;
  double mlo, mhi, dlogm, nsat, m1, mbar, slogm=0.7, ptot, phalo, pn, nhalo, mass, p[401];

  m1 = pow(10.0,12.85);

  mbar = pow(10.0,13.85);
  mlo = mbar/30;
  mhi = mbar*30;
  nmass = 1000;
  dlogm = log(mhi/mlo)/nmass;

  nmax = 30;
  fprintf(stderr,"%d\n",nmax);

  ptot = 0;
  for(n=0;n<=nmax;++n)
    {
      p[n] = 0;
      for(i=0;i<nmass;++i)
	{
	  mass = exp((i+0.5)*dlogm )*mlo;
	  nhalo = halo_abundance(mass)*mass;
	  phalo = 1/sqrt(2*PI)/slogm*exp(-log(mass/mbar)*log(mass/mbar)/(2*slogm*slogm));
	  nsat = mass/m1;
	  pn = poisson_prob(n,nsat);
	  p[n] += dlogm*nhalo*phalo*pn;
	  ptot += dlogm*nhalo*phalo*pn;
	}
    }
  for(n=0;n<=nmax;++n)
    printf("%d %e %e\n",n,p[n]/ptot,poisson_prob(n,mbar/m1));
}

/* Poisson probability of n given n_average
 */
double poisson_prob(int n, double nave)
{
  int i;
  double fac=1;

  if(n>0)
    for(i=1;i<=n;++i)
      fac*=nave/i;

  return((float)(fac*exp(-nave)));
}

