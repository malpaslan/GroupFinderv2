#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

#define HALO_MAX 1.0E+16

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

/* Local functions
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
float func_match_nhost(float mass);

float density2halo(float galaxy_density)
{
  g7_ngal = galaxy_density;
  return exp(zbrent(func_match_ngal,log(1.0E+8),log(HALO_MAX),1.0E-5));
  
}
float density2host_halo(float galaxy_density)
{
  g7_ngal = galaxy_density;
  return exp(zbrent(func_match_nhost,log(1.0E+8),log(HALO_MAX),1.0E-5));
  
}

float func_match_ngal(float mass)
{
  static int flag=1, n=100, prev_cosmo=-1, call_count=1;
  static float *mh, *ms, *mx, *nh, mlo, mhi, dlogm, mmax;
  int i;
  float a, maglo, maghi, dmag, m, n1, n2;

  
  if(flag)
    {
      flag = 0;
      mh = vector(1,n);
      nh = vector(1,n);
      mx = vector(1,n);

      mlo = 1.0E+8;
      mhi = HALO_MAX;
      dlogm = log(mhi/mlo)/n;

      for(i=1;i<=n;++i)
	{
	  mh[i] = exp((i-0.5)*dlogm)*mlo;
	  n1 = qromo(halo_abundance2,log(mh[i]),log(HALO_MAX),midpnt);
	  if(mh[i]<HALO_MAX/1.0)
	    n2 = qromo(subhalo_abundance2,log(mh[i]),log(HALO_MAX/1.0),midpnt);
	  else
	    n2 = 0;
	  if(n2<0)n2=0;
	  nh[i] = log(n1+n2);
	  //nh[i] = log(qromo(func_nhalo,log(mh[i]),log(HALO_MAX),midpnt));
	  //printf("SHAM %e %e %e %e %e\n",(mh[i]),exp(nh[i]),n1,n2,subhalo_abundance2(log(mh[i])));
	  mh[i] = log(mh[i]);
	  fflush(stdout);
	}
      spline(mh,nh,n,1.0E+30,1.0E+30,mx);
    }
  splint(mh,nh,mx,n,mass,&a);
  return exp(a)-g7_ngal;
}

float func_match_nhost(float mass)
{
  static int flag=1, n=100, prev_cosmo=-1, call_count=1;
  static float *mh, *ms, *mx, *nh, mlo, mhi, dlogm, mmax;
  int i;
  float a, maglo, maghi, dmag, m, n1, n2;

  
  if(flag)
    {
      flag = 0;
      mh = vector(1,n);
      nh = vector(1,n);
      mx = vector(1,n);

      mlo = 1.0E+8;
      mhi = HALO_MAX;
      dlogm = log(mhi/mlo)/n;

      for(i=1;i<=n;++i)
	{
	  mh[i] = exp((i-0.5)*dlogm)*mlo;
	  n1 = qromo(halo_abundance2,log(mh[i]),log(HALO_MAX),midpnt);
	  nh[i] = log(n1);
	  mh[i] = log(mh[i]);
	  fflush(stdout);
	}
      spline(mh,nh,n,1.0E+30,1.0E+30,mx);
    }
  splint(mh,nh,mx,n,mass,&a);
  return exp(a)-g7_ngal;
}


// recall that the input mmin is already in log units
float func_match_nhalo(float mmin)
{
  return qromo(func_nhalo,mmin,log(HALO_MAX),midpnt) - g7_ngal;
}

float func_nhalo(float m)
{
  m = exp(m);
  //printf("%e %e\n",m,subhalo_abundance(m));
  return (halo_abundance(m)+subhalo_abundance(m))*m;
}

/* integrate over the parent halo mass function
 * to get the density of subhalos of mass m
 */
float subhalo_abundance(float m)
{  
  g7_msub = m;
  return qromo(func_subhalo1,log(m),log(HALO_MAX),midpnt);
}

float subhalo_abundance2(float m)
{  
  g7_msub = exp(m);
  //printf("%e\n",g7_msub);
  return qromo(func_subhalo1,(m),log(HALO_MAX),midpnt)*g7_msub;
}


/* NB the 1.5 factor is to fit the results from the z=0 200 box
 */
float func_subhalo1(float mhost)
{
  double x;
  mhost = exp(mhost);
  x= pow(g7_msub/mhost,-0.7)*exp(-9.9*pow(g7_msub/mhost,2.5))*0.3*halo_abundance(mhost)*mhost/g7_msub;
  //printf("%e %e\n",mhost,x)
  return x;
  //  return pow(g7_msub/mhost,-0.8)*exp(-g7_msub/mhost*1.25)*0.2*halo_abundance(mhost)*mhost/g7_msub;
}
 
float halo_abundance2(float m)
{
  m=exp(m);
  return halo_abundance(m)*m;
}

float halo_abundance(float m)
{
  int i;
  FILE *fp;
  float a;
  static int n=0;
  static float *x, *y, *z;

  if(!n)
    {
      fp = openfile("censat_dndm.dat");
      //fp = openfile("wmap1.massfunc");
      //fp = openfile("s8_0.7.massfunc");
      n = filesize(fp);
      x = vector(1,n);
      y = vector(1,n);
      z = vector(1,n);
      for(i=1;i<=n;++i)
	{
	  fscanf(fp,"%f %f",&x[i],&y[i]);
	  x[i] = log(x[i]);
	  y[i] = log(y[i]);
	}
      spline(x,y,n,1.0E+30,1.0E+30,z);
      fclose(fp);
    }
  splint(x,y,z,n,log(m),&a);
  return exp(a);
}

float func_ngal(float mag)
{
  float phi1=0.0156;
  float phi2=0.0062;
  float mstar=-20.04;
  float a1=-0.17;
  float a2=-1.52;

  // temp! replacing with blanton 03 for a second...
  return 0.4*2.30258*exp(-pow(10.0,-0.4*(mag +20.44)))*(1.49e-2*pow(10.0,-0.4*(mag + 20.44)*(-1.05+1)));

  // blanton's 05 lowL LF
  return 0.4*2.30258*exp(-pow(10.0,-0.4*(mag - mstar)))*
    (phi1*pow(10.0,-0.4*(mag - mstar)*(a1+1)) + phi2*pow(10.0,-0.4*(mag-mstar)*(a2+1)));

}
