// gcc -o popHalo populate_halo.c *.o -L/home/users/ma5046/libC_main/ -lC_main -lm
// ~tinker/exec/wp_covar 0.1 20 10 250 0 250 1 ../misc_work/bolshoi_halopop_10.9.dat a 0 1 1 1 5


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <sys/time.h>
#include "nrutil.h"
#include "header.h"

// Definitions

#define OMEGA_M 0.25
#define RHO_CRIT 2.775E+11
#define DELTA_HALO 200
#define OMEGA_TEMP 0.3
#define CVIR_FAC 1.0
#define VBIAS 1.0
#define VBIAS_C 0.0
#define REDSHIFT 0
#define MCMC 0
#define BOX_SIZE 250.0

// Internal functions go here:

char *concat(const char *s1, const char *s2);
int poisson_deviate(double nave);
double NFW_central_velocity(double mass, double v[], double mag);
double NFW_density(double r, double rs, double ps);
double NFW_velocity(double mass, double v[], double mag);
double NFW_position(double mass, double x[]);

// This code populates halos from an external file (output from a simulation) with mock galaxies, assuming some HOD. Normally one would parametrize this, but this version will read in a HOD file.
// Input argument to be one of the following:
// 9.5  9.8 10.1 10.3 10.6 10.9 11.2 11.4 11.7 12.0

int main(int argc, char **argv){

	FILE *fp, *fpoutQ, *fpoutNQ;
	int haloID, i, j, k, l, n, count, counter ,nhalo, HODQdim, HODNQdim, numhalo, n1;
	float *mass, **poshalo, **velhalo, ncen, nsat, **HODQcen, **HODQsat, **HODNQcen, **HODNQsat;
	float mMinCQ, mMinSQ, mMaxCQ, mMaxSQ, nMaxSQ, nMaxCQ, nMinSQ, nMinCQ;
  float mMinCNQ, mMinSNQ, mMaxCNQ, mMaxSNQ, nMaxSNQ, nMaxCNQ, nMinSNQ, nMinCNQ;
  double r, mag, posgal[3], velgal[3], tempmass;
  float x, x0, x1, y, y0, y1;
	char *ff, ffoutQ[10000],ffoutNQ[10000],tempff[1000];
	char string[1000];
  long IDUM = -445;

	// Read in the halo file.
  count = 0;
	ff = "/home/users/ma5046/misc_work/bolshoi_redux.dat";
  
 	fp = fopen(ff, "r");
	if(!(fp=fopen(ff,"r"))){
      printf("ERROR opening [%s]\n",ff);
      exit(0);
    }

    // Measure length of file.

    while(fgets(string,1000,fp)){
  		count++;
	}
	rewind(fp);

	nhalo = count-1; // First line of file contains column headers.
	count = 0;
	printf("%d halos.\n",nhalo);

	// Load halo information into arrays / matrices.

	mass = vector(1, nhalo);
	//ncen = vector(1, nhalo);
	//nsat = vector(1, nhalo);
	poshalo = matrix(1,nhalo,1,3);
	velhalo = matrix(1,nhalo,1,3);

	fgets(string,1000,fp);

  // Read mass into a temporary variable and store to the mass array after logging.

	for(i = 1; i <= nhalo; ++i){
		fscanf(fp, "%f %f %f %f %f %f %f", &poshalo[i][1], &poshalo[i][2], &poshalo[i][3],&velhalo[i][1], &velhalo[i][2], &velhalo[i][3],&mass[i]);
	}	
	
	fclose(fp);
  
 	// Data read in! Now loop through halos and compute number of satellites and centrals.
	
	// Read in the HOD specified in the command. This is probably a hackneyed way of doing this in C, but I am also very bad at C.
	
  // Start with quenched galaxies.
	// Central HOD.
  char fftab[1000];

  sprintf(tempff,"cenQTable%s.dat",argv[1]);
  sprintf(fftab,"/home/users/ma5046/misc_work/halotables/%s",tempff);
	printf("%s\n",fftab);
  fp = fopen(fftab, "r");
	while(fgets(string,1000,fp)){
		count++;
	}
	rewind(fp);
	
	HODQdim = count;
  printf("%d\n",count);
	count = 0;
	HODQcen = matrix(1,HODQdim,1,2);
				
	for(i = 1; i <= HODQdim; ++i){
		fscanf(fp,"%f %f", &HODQcen[i][1], &HODQcen[i][2]);	
	}
	fclose(fp);
	mMinCQ = HODQcen[1][1];
  nMinCQ = HODQcen[1][2];
	mMaxCQ = HODQcen[HODQdim][1];
  nMaxCQ = HODQcen[HODQdim][2];

	// Satellite HOD.

	sprintf(tempff,"satQTable%s.dat",argv[1]);
  sprintf(fftab, "/home/users/ma5046/misc_work/halotables/%s",tempff);
  printf("%s\n",fftab);
  fp = fopen(fftab, "r");
  while(fgets(string,1000,fp)){
    count++;
  }
  rewind(fp);

  HODQdim = count;
  HODQsat = matrix(1,HODQdim,1,2);
	count = 0;
  for(i = 1; i <= HODQdim; ++i){
    fscanf(fp,"%f %f", &HODQsat[i][1], &HODQsat[i][2]);
  }
  
	fclose(fp);
	mMinSQ = HODQsat[1][1];
  nMinSQ = HODQsat[1][2];
	mMaxSQ = HODQsat[HODQdim][1];
  nMaxSQ = HODQsat[HODQdim][2];

	printf("Quenched: Stellar mass = %s\nCentral min halo mass = %f; max halo mass = %f\nSatellite min halo mass = %f; max halo mass = %f\n\n",argv[1],mMinCQ,mMaxCQ,mMinSQ,mMaxSQ);

  // Time for unquenched galaxies.
  // Central HOD.
  sprintf(tempff,"cenNQTable%s.dat",argv[1]);	
  sprintf(fftab,"/home/users/ma5046/misc_work/halotables/%s",tempff);
	printf("%s\n",fftab);
  fp = fopen(fftab, "r");
	while(fgets(string,1000,fp)){
		count++;
	}
	rewind(fp);
	
	HODNQdim = count;
  printf("%d\n",count);
	count = 0;
	HODNQcen = matrix(1,HODNQdim,1,2);
				
	for(i = 1; i <= HODNQdim; ++i){
		fscanf(fp,"%f %f", &HODNQcen[i][1], &HODNQcen[i][2]);	
	}
	fclose(fp);
	mMinCNQ = HODNQcen[1][1];
  nMinCNQ = HODNQcen[1][2];
	mMaxCNQ = HODNQcen[HODNQdim][1];
  nMaxCNQ = HODNQcen[HODNQdim][2];

	// Satellite HOD.

	sprintf(tempff,"satNQTable%s.dat",argv[1]);
  sprintf(fftab, "/home/users/ma5046/misc_work/halotables/%s",tempff);
  printf("%s\n",fftab);
  fp = fopen(fftab, "r");
  while(fgets(string,1000,fp)){
    count++;
  }
  rewind(fp);

  HODNQdim = count;
  HODNQsat = matrix(1,HODNQdim,1,2);
	count = 0;
  for(i = 1; i <= HODNQdim; ++i){
    fscanf(fp,"%f %f", &HODNQsat[i][1], &HODNQsat[i][2]);
  }
  
	fclose(fp);
	mMinSNQ = HODNQsat[1][1];
  nMinSNQ = HODNQsat[1][2];
	mMaxSNQ = HODNQsat[HODNQdim][1];
  nMaxSNQ = HODNQsat[HODNQdim][2];

	printf("Not quenched: Stellar mass = %s\nCentral min halo mass = %f; max halo mass = %f\nSatellite min halo mass = %f; max halo mass = %f\n\n",argv[1],mMinCNQ,mMaxCNQ,mMinSNQ,mMaxSNQ);

	// Time to populate the halo.

	printf("Data read in! Time to populate the halo...\n");
 
  // Set up output file.
  sprintf(ffoutQ,"/home/users/ma5046/misc_work/bolshoi_Qhalopop_%s.dat",argv[1]);
  sprintf(ffoutNQ,"/home/users/ma5046/misc_work/bolshoi_NQhalopop_%s.dat",argv[1]);
  
  printf("%s\n",ffoutQ);
  printf("%s\n",ffoutNQ);
  fpoutQ = fopen(ffoutQ, "w");
  fpoutNQ = fopen(ffoutNQ,"w");

  // Loop through every halo. If its mass is lower than the minimum mass, it has no central/satellite. If this isn't the case, linearly interpolate the HOD values to get number of centrals and satellites.

  haloID = 0;

// Start with the quenched galaxies.

	for(i = 1; i <= nhalo; ++i){

    haloID++;
    
	  tempmass = (double) pow(10.0,mass[i]);
    
	// Centrals.
		if(mass[i] < mMinCQ | mass[i] > mMaxCQ){
			goto SATELLITESQ;
		}
	  	
		// Identify what rows to interpolate between in the HOD.
		k = 1;
		while(HODQcen[k][1] < mass[i]){
      k++;
    }
    // Interpolation will be borked if k is the final row of the array. Account for this.
 		
		// k is the row.
    if(HODQcen[k][1] != mMinCQ){   
      x1 = HODQcen[k][1];
      y1 = HODQcen[k][2];
      x0 = HODQcen[k-1][1];
      y0 = HODQcen[k-1][2];
    }
    else{
      x1 = HODQcen[k+1][1];
      y1 = HODQcen[k+1][2];
      x0 = mMinCQ;
      y0 = nMinCQ;
    }
    x = mass[i];
    y = (y0*(x1-x) + y1*(x-x0))/(x1-x0);
    
    // Reset k for next halo.
    
    k = 1;
  
    ncen = pow(10.0,y);
    
    //ncen = y;
    
    if(ncen > 1.0) ncen = 1.0;
    
    // Draw a random number. If greater than ncen, go to satellites. (So if ncen = 1, always generate a central).
    
    if(drand48() > ncen) goto SATELLITESQ;
    
    // insert NFW position and velocity stuff here.
    
    NFW_central_velocity(tempmass,velgal,mag);
    //printf("cen %d\n",i);
    for(j = 0; j <= 2; ++j){
      // Add NFW velocity of galaxy to halo's velocity.
      velhalo[i][j] += velgal[j];
    }

    fprintf(fpoutQ,"%f %f %f %f %f %f %f %d\n",poshalo[i][0],poshalo[i][1],poshalo[i][2],velhalo[i][0],velhalo[i][1],velhalo[i][2],mass[i], 1);

    //printf("ncen: %f %d %e\n",y,numhalo,ncen);
   
	  // Satellites.
	  SATELLITESQ:
    if(mass[i] < mMinSQ | mass[i] > mMaxSQ){
      continue;
    }
    

    k = 1;
    while(HODQsat[k][1] < mass[i]){
      ++k;
    }
//printf("sat %f %f %f\n",HODsat[k][1], mass[i], HODsat[k+1][1]);	
    if(HODQsat[k][1] != mMinSQ){   
      x1 = HODQsat[k][1];
      y1 = HODQsat[k][2];
      x0 = HODQsat[k-1][1];
      y0 = HODQsat[k-1][2];
    }
    else{
      x1 = HODQsat[k+1][1];
      y1 = HODQsat[k+1][2];
      x0 = mMinSQ;
      y0 = nMinSQ;
    }
    x = mass[i];
    y = (y0*(x1-x) + y1*(x-x0))/(x1-x0);
    nsat = pow(10.0,y);
    //nsat = y;
    k = 1;
    
    // Create a poisson deviate for nsat.
		if(nsat > 250){
      n1 = gasdev(&IDUM)*sqrt(nsat)+nsat;
    }
    else
      n1 = poisson_deviate(nsat);
    
    // Loop through each satellite and give it a position and velocity.

    for(j = 1; j <= n1; ++j){
      //printf("sat %d %d %lf\n", j, i, tempmass);
     
      r = NFW_position(tempmass,posgal);
      NFW_velocity(tempmass,velgal,mag);
      for(k = 0; k <= 2; ++k){
        // Modify galaxy position and velocity. Careful for box edges.

        posgal[k] += poshalo[i][k];
        if(posgal[k] < 0) posgal[k] += BOX_SIZE;
        else if(posgal[k] > BOX_SIZE) posgal[k] -= BOX_SIZE;
        velgal[k] += velhalo[i][k];
      }
      
      fprintf(fpoutQ,"%f %f %f %f %f %f %f %d\n",posgal[0],posgal[1],posgal[2],velgal[0],velgal[1],velgal[2],mass[i], 0);
    }
    //printf("nsat %f\n:", nsat);
  }
  fclose(fpoutQ);

  // Unquenched galaxies now!
  
  for(i = 1; i <= nhalo; ++i){

    haloID++;
    
	  tempmass = (double) pow(10.0,mass[i]);
    
	  // Centrals.
		if(mass[i] < mMinCNQ | mass[i] > mMaxCNQ){
			goto SATELLITESNQ;
		}
	  	
		// Identify what rows to interpolate between in the HOD.
		k = 1;
		while(HODNQcen[k][1] < mass[i]){
      k++;
    }
    // Interpolation will be borked if k is the final row of the array. Account for this.
 		
		// k is the row.
    if(HODNQcen[k][1] != mMinCNQ){   
      x1 = HODNQcen[k][1];
      y1 = HODNQcen[k][2];
      x0 = HODNQcen[k-1][1];
      y0 = HODNQcen[k-1][2];
    }
    else{
      x1 = HODNQcen[k+1][1];
      y1 = HODNQcen[k+1][2];
      x0 = mMinCNQ;
      y0 = nMinCNQ;
    }
    x = mass[i];
    y = (y0*(x1-x) + y1*(x-x0))/(x1-x0);
    
    // Reset k for next halo.
    
    k = 1;
  
    ncen = pow(10.0,y);
    
    //ncen = y;
    
    if(ncen > 1.0) ncen = 1.0;
    
    // Draw a random number. If greater than ncen, go to satellites. (So if ncen = 1, always generate a central).
    
    if(drand48() > ncen) goto SATELLITESNQ;
    
    // insert NFW position and velocity stuff here.
    
    NFW_central_velocity(tempmass,velgal,mag);
    //printf("cen %d\n",i);
    for(j = 0; j <= 2; ++j){
      // Add NFW velocity of galaxy to halo's velocity.
      velhalo[i][j] += velgal[j];
    }

    fprintf(fpoutNQ,"%f %f %f %f %f %f %f %d\n",poshalo[i][0],poshalo[i][1],poshalo[i][2],velhalo[i][0],velhalo[i][1],velhalo[i][2],mass[i], 1);

    //printf("ncen: %f %d %e\n",y,numhalo,ncen);
   
	  // Satellites.
	  SATELLITESNQ:
    if(mass[i] < mMinSNQ | mass[i] > mMaxSNQ){
      continue;
    }
    

    k = 1;
    while(HODNQsat[k][1] < mass[i]){
      ++k;
    }
//printf("sat %f %f %f\n",HODsat[k][1], mass[i], HODsat[k+1][1]);	
    if(HODNQsat[k][1] != mMinSNQ){   
      x1 = HODNQsat[k][1];
      y1 = HODNQsat[k][2];
      x0 = HODNQsat[k-1][1];
      y0 = HODNQsat[k-1][2];
    }
    else{
      x1 = HODNQsat[k+1][1];
      y1 = HODNQsat[k+1][2];
      x0 = mMinSNQ;
      y0 = nMinSNQ;
    }
    x = mass[i];
    y = (y0*(x1-x) + y1*(x-x0))/(x1-x0);
    nsat = pow(10.0,y);
    //nsat = y;
    k = 1;
    
    // Create a poisson deviate for nsat.
		if(nsat > 250){
      n1 = gasdev(&IDUM)*sqrt(nsat)+nsat;
    }
    else
      n1 = poisson_deviate(nsat);
    
    // Loop through each satellite and give it a position and velocity.

    for(j = 1; j <= n1; ++j){
      //printf("sat %d %d %lf\n", j, i, tempmass);
     
      r = NFW_position(tempmass,posgal);
      NFW_velocity(tempmass,velgal,mag);
      for(k = 0; k <= 2; ++k){
        // Modify galaxy position and velocity. Careful for box edges.

        posgal[k] += poshalo[i][k];
        if(posgal[k] < 0) posgal[k] += BOX_SIZE;
        else if(posgal[k] > BOX_SIZE) posgal[k] -= BOX_SIZE;
        velgal[k] += velhalo[i][k];
      }
      
      fprintf(fpoutNQ,"%f %f %f %f %f %f %f %d\n",posgal[0],posgal[1],posgal[2],velgal[0],velgal[1],velgal[2],mass[i], 0);
    }
    //printf("nsat %f\n:", nsat);
  }

  fclose(fpoutNQ);
	exit(0);
}

/* Poisson probability of n given n_average
 *  */
//double poisson_prob(int n, double nave){
//  int i;
//  double fac=1;
//
//  if(n>0)
//    for(i=1;i<=n;++i)
//      fac*=nave/i;
//
//  return((float)(fac*exp(-nave)));
//}

/* Randomly generates a position away from the origin with 
 *  * a probability given by the NFW profile for a halo of the input
 *   * mass (and including the CVIR_FAC)
 *    */
double NFW_position(double mass, double x[]){
  double r,pr,max_p,costheta,sintheta,phi1,signs,rvir,rs,cvir, mfac=1;

  cvir=halo_concentration(mass);
  rvir=pow(3*mass*mfac/(4*DELTA_HALO*PI*RHO_CRIT*OMEGA_M),1.0/3.0);
  rs=rvir/cvir;
  max_p=NFW_density(rs,rs,1.0)*rs*rs*4.0*PI;

  for(;;) {
    r=drand48()*rvir;
    pr=NFW_density(r,rs,1.0)*r*r*4.0*PI/max_p;
                              
    if(drand48()<=pr){
      costheta=2.*(drand48()-.5);
      sintheta=sqrt(1.-costheta*costheta);
      signs=2.*(drand48()-.5);
      costheta=signs*costheta/fabs(signs);
      phi1=2.0*PI*drand48();
                                                  	
      x[0]=r*sintheta*cos(phi1);
      x[1]=r*sintheta*sin(phi1);
      x[2]=r*costheta;
      return r;
    }
  }
}

/* This is the NFW density profile
 *  */
double NFW_density(double r, double rs, double ps){
    return(ps*rs/(r*(1+r/rs)*(1+r/rs)));
}

/* This sets the velocity to be isotropic Gaussian.
 *  */
double NFW_velocity(double mass, double v[], double mag){
  static long IDUM2=-455;
  static double fac = -1;
  double sigv,vbias=1,mfac=1;
  int i;

  if(MCMC>=3)
    mfac = 0.30711/OMEGA_TEMP;
  if(fac<0)
    fac=sqrt(4.499E-48)*pow(4*DELTA_HALO*PI*OMEGA_M*RHO_CRIT/3,1.0/6.0)*3.09E19*sqrt(1+REDSHIFT);
  sigv=fac*pow(mass*mfac,1.0/3.0)/sqrt(2.0);
  for(i=0;i<3;++i)
    v[i]=gasdev(&IDUM2)*sigv*VBIAS;
  return(0);
}

/* This sets the velocity to be isotropic Gaussian.
 *  */
double NFW_central_velocity(double mass, double v[], double mag){
  static long IDUM2=-455;
  static double fac = -1;
  double sigv,vbias=1;
  int i;

  if(fac<0)
    fac=sqrt(4.499E-48)*pow(4*DELTA_HALO*PI*OMEGA_M*RHO_CRIT/3,1.0/6.0)*3.09E19;
  sigv=fac*pow(mass,1.0/3.0)/sqrt(2.0);
  for(i=0;i<3;++i)
    v[i]=gasdev(&IDUM2)*sigv*VBIAS_C;
  return(0);
}

char* concat(const char *s1, const char *s2){  
    char *result = malloc(strlen(s1)+strlen(s2)+1);//+1 for the null-terminator
    strcpy(result, s1);
    strcat(result, s2);
    return result;
}

/* Generate a random integer based on a Poisson distribution 
 *  * with mean given as input. Code by JLT.
 *   */
int poisson_deviate(double nave){
  static int flag=0;
  double p,pp;
  int n;

  p=0;
  pp=1;

  while(p<pp){
    if(nave<1)
      n=(int)(drand48()*20);
    else
      n=(int)(drand48()*30*nave);
    p=poisson_prob(n,nave);
    pp=drand48();
    }
    return(n);
}
