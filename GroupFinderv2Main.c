// GroupFinderv2 by Mehmet Alpaslan.
// Adapted from Jeremy Tinker's group finding algorithm (Tinker et al. 2011).

// Principal changes include modifications to improve run time; chi-squared optimization; and sensitivity to galaxy-halo mass relation scatter. Other tweaks also.

// Redshifts should be kept in velocity space for computations (for now).

// Initialization //

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

// Definitions

#define omegaM 0.25
#define pi 3.141592741
#define rhoCrit 2.775E+11
#define dHalo 200
#define speedOfLight 3.0E+5
#define cOnH0 2997.92
#define bigG 4.304E-9 /* BIG G in units of (km/s)^2*Mpc/M_sol */
#define g0 (1.0/sqrt(2.0*3.14159))
#define root2 1.41421
#define q0 2.0
#define q1 -1.0
#define qz0 0.1
#define third (1.0/3.0)
#define ang (pi/180.0)
#define rt2Pi 2.50663
#define czMin 7000.0
#define czBuf 0
//#define REDSHIFT (12000.0/SPEED_OF_LIGHT)// -18
//#define MAGNITUDE -18.0
//#define REDSHIFT (19200.0/SPEED_OF_LIGHT)// -19
//#define MAGNITUDE -19.0
// #define REDSHIFT (31000.0/SPEED_OF_LIGHT)// -20
// #define MAGNITUDE -20.0

// Imported functions from numerical recipes 

float qromo(float (*func)(float), float a, float b,
	     float (*choose)(float(*)(float), float, float, int));
float midpnt(float (*func)(float), float a, float b, int n);
void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);
void sort2(int n, float arr[], int id[]);
void sort3(int n, float arr[], float brr[]);
float qtrap(float (*func)(float), float a, float b);
float *vector(long nl, long nh);
int *ivector(long nl, long nh);

// External functions go here

float fiber_corrected_galaxy_property(int indx, float *mag_r, float *mag_g, float *prop, int ngal, int *flag);
float density2halo(float galaxy_density);
float density2host_halo(float galaxy_density);

// Local functions go here

float distance_redshift(float z);
float angular_separation(float a1, float d1, float a2, float d2);
float find_satellites(int i, float *ra, float *dec, float *redshift, float *mag_r, float theta_max, float x1, int *group_member, int *indx, int ngal, float radius, float mass, int igrp, float *luminosity, float *nsat_cur, int i1, float *prob_total);
float radial_probability(float mass, float dr, float rad, float ang_rad);

// Global variables go here

float GALAXY_DENSITY, MSTARLIM, MAGNITUDE, MAXREDSHIFT, MINREDSHIFT;

// End initialization //

// Here we go!

int main(int argc, char **argv){

	int i, j, k, igrp, ngal, nsample, count;
	int *indx, *collision, *ka;
	int *group_member, *group_index, *group_center, *temp_group;
	int nsat_tot, ngrp, xngrp, niter, niter_max, ngrp_temp;
	float *ra, *dec, *redshift, *mag_g, *mag_r, *v_max, *m_stellar, *Hdelta, *Dn4000, *luminosity, *Rexp, *sSFR, *sersic, *velDisp, *sNr, *petroRad;
	float *mass, *rad, *angRad, *sigma, *prob_total, *nsat_indi;
	float *group_luminosity;
	float x1, x2, *tempArray;
	float volume;
	double ndens_gal = 0;
	char string[1000];
	char *ff;
	FILE *fp;

	count = 0;
	niter_max = 10;

	MAXREDSHIFT = 19200/speedOfLight;
	MINREDSHIFT = czMin / speedOfLight;
	MAGNITUDE = -19;
	MSTARLIM = pow(10.0,9.8);

	volume = 4./3.*pi*(pow(distance_redshift(MAXREDSHIFT),3.0)-pow(distance_redshift(MINREDSHIFT),3.0))*0.177;

	// Start by reading in data. At first only read in bare minimum to establish a volume limited sample; then read in everything else.

	fprintf(stderr,"\n** Reading in data and defining sample. **\n");

	// Import VAGC main catalogue (RA, Dec, z). 
	// Measure length of this catalogue; this is the total number of galaxies. Apply this length to the arrays containing redshift, magnitude, and mass.

	ff = "/Users/mehmet/Dropbox/nyucosmo/vagc_files_groups/lss.dr72bright34.dat";

	fp = fopen(ff,"r");
	if(!(fp=fopen(ff,"r"))){
      fprintf(stderr,"ERROR opening [%s]\n",ff);
      exit(0);
    }

    // This measures the length of the file.

	while(fgets(string,1000,fp)){
  		count++;
	}
	rewind(fp);
	// Assign count to ngal; this is the number of galaxies.

	ngal = count;

	// Each variable that holds a different galaxy property can now be made into an array going from 1 to ngal.

	ra = vector(1,ngal);
  	dec = vector(1,ngal);
  	redshift = vector(1,ngal);
  	mag_r = vector(1, ngal);
  	mag_g = vector(1, ngal);
  	m_stellar = vector(1, ngal);
  	indx = ivector(1, ngal);

  	for(i = 1; i <= ngal; ++i){
  		fscanf(fp,"%*d %*d %*d %f %f %f %*f %*f",&ra[i],&dec[i],&redshift[i]);
  		ra[i] *= pi/180;
  		dec[i] *= pi/180;
  		indx[i] = i;
  		redshift[i] /= 3E5;
  		fgets(string,1000,fp);
  	}

  	fclose(fp);
  	fprintf(stderr, "Astrometry read in.\n");

  	// Import photometry (mag_r, mag_g).

  	ff = "/Users/mehmet/Dropbox/nyucosmo/vagc_files_groups/photo_evolve.dr72bright34.dat";

	fp = fopen(ff,"r");
	if(!(fp=fopen(ff,"r"))){
    	fprintf(stderr,"ERROR opening [%s]\n",ff);
    	exit(0);
    }

    for(i = 1; i <= ngal; ++i){
    	fscanf(fp, "%*d %*f %f %f %*f %*f %*f %*f", &mag_g[i],&mag_r[i]);
    }  	

    fclose(fp);
    fprintf(stderr, "Photometry read in.\n");

    // Import stellar mass (m_stellar).

    ff = "/Users/mehmet/Dropbox/nyucosmo/vagc_files_groups/mstellar.dr72bright34.dat";

	fp = fopen(ff,"r");
	if(!(fp=fopen(ff,"r"))){
    	fprintf(stderr,"ERROR opening [%s]\n",ff);
    	exit(0);
    }

    for(i = 1; i <= ngal; ++i){
    	fscanf(fp, "%*d %f", &m_stellar[i]);
    }

    fclose(fp);
    fprintf(stderr,"Stellar masses read in.\n");

    // Now have enough data to make a volume limited sample.
    // Start by identifying how many galaxies are in this sample. Assign this length to all remaining data arrays.

    count = 0;
    for(i = 1; i <= ngal; ++i){
    	if(mag_r[i] <= MAGNITUDE && redshift[i] <= MAXREDSHIFT && redshift[i] >= MINREDSHIFT && m_stellar[i] >= MSTARLIM){
    		count++;
    	}
    }

    nsample = count;

    v_max = vector(1, nsample);
  	luminosity = vector(1, nsample);
  	Hdelta = vector(1, nsample);
  	ka = ivector(1, nsample);
  	Rexp = vector(1, nsample);
  	Dn4000 = vector(1, nsample);
  	sSFR = vector(1, nsample);
  	sersic = vector(1, nsample);
  	collision = ivector(1, nsample);
  	velDisp = vector(1, nsample);
  	sNr = vector(1, nsample);
  	petroRad = vector(1, nsample);

  	mass = vector(1, nsample);
  	rad = vector(1, nsample);
  	angRad = vector(1, nsample);
  	sigma = vector(1, nsample);
  	prob_total = vector(1, nsample);

  	tempArray = vector(1, nsample);
  	temp_group = ivector(1, nsample);

  	group_luminosity = vector(1, nsample);
  	nsat_indi = vector(1, nsample);
  	group_member = ivector(1, nsample);
  	group_index = ivector(1, nsample);
  	group_center = ivector(1, nsample);
  	
    // Import Vmax (v_max).

    ff = "/Users/mehmet/Dropbox/nyucosmo/vagc_files_groups/vmax_evolve.dr72bright34.dat";

    fp = fopen(ff,"r");
	if(!(fp=fopen(ff,"r"))){
    	fprintf(stderr,"ERROR opening [%s]\n",ff);
    	exit(0);
    }

    j = 1;
    for(i = 1; i <= ngal; ++i){
    	if(mag_r[i] <= MAGNITUDE && redshift[i] <= MAXREDSHIFT && redshift[i] >= MINREDSHIFT && m_stellar[i] >= MSTARLIM){
    		fscanf(fp, "%*d %f %*f %*f", &v_max[j]);
    		++j;
    	}
    }
    fclose(fp);

    if(j-1 != nsample){
    	printf("Error! v_max array size not the same as sample size: %d %d\n", j-1, nsample);
    	exit(0);
    }

    fprintf(stderr,"Vmax read in.\n");

    // Import Hdelta (Hdelta).

    ff = "/Users/mehmet/Dropbox/nyucosmo/vagc_files_groups/Hdelta_A.dr72bright34.dat";

	fp = fopen(ff,"r");
	if(!(fp=fopen(ff,"r"))){
    	fprintf(stderr,"ERROR opening [%s]\n",ff);
    	exit(0);
    }    

    j = 1;
    for(i = 1; i <= ngal; ++i){
    	if(mag_r[i] <= MAGNITUDE && redshift[i] <= MAXREDSHIFT && redshift[i] >= MINREDSHIFT && m_stellar[i] >= MSTARLIM){
    		fscanf(fp, "%f %*f %*d %*d %*f %*f %*f %*f %*f", &Hdelta[j]);
    		++j;
    	}
    }
    fclose(fp);

    if(j-1 != nsample){
    	printf("Error! Hdelta array size not the same as sample size: %d %d\n", j-1, nsample);
    	exit(0);
    }

    fprintf(stderr, "Hdelta read in.\n");

    // Import KA parameters (ka).

    ff = "/Users/mehmet/Dropbox/nyucosmo/vagc_files_groups/ka.dr72bright34.dat";

    fp = fopen(ff,"r");
	if(!(fp=fopen(ff,"r"))){
    	fprintf(stderr,"ERROR opening [%s]\n",ff);
    	exit(0);
    }    

    j = 1;
    for(i = 1; i <= ngal; ++i){
    	if(mag_r[i] <= MAGNITUDE && redshift[i] <= MAXREDSHIFT && redshift[i] >= MINREDSHIFT && m_stellar[i] >= MSTARLIM){
    		fscanf(fp, "%f %f %*f", &x1, &x2);
    		ka[j] = 0;
    		if(x1 > 0.2 && x2 < pow(10.0,x1) - 1)
    			ka[j] = 1;
    		++j;
    	}
    }
    fclose(fp);

    if(j-1 != nsample){
    	printf("Error! ka array size not the same as sample size: %d %d\n", j-1, nsample);
    	exit(0);
    }

    fprintf(stderr, "KA parameters read in.\n");

    // Import galaxy radius (R_exp).

	ff = "/Users/mehmet/Dropbox/nyucosmo/vagc_files_groups/galrad.dr72bright34.dat";

    fp = fopen(ff,"r");
	if(!(fp=fopen(ff,"r"))){
    	fprintf(stderr,"ERROR opening [%s]\n",ff);
    	exit(0);
    }        

    j = 1;
    for(i = 1; i <= ngal; ++i){
    	if(mag_r[i] <= MAGNITUDE && redshift[i] <= MAXREDSHIFT && redshift[i] >= MINREDSHIFT && m_stellar[i] >= MSTARLIM){
    		fscanf(fp, "%f %*f %*f %*f %*f", &Rexp[j]);
    		++j;
    	}
    }
	fclose(fp);

	if(j-1 != nsample){
    	printf("Error! R_exp array size not the same as sample size: %d %d\n", j-1, nsample);
    	exit(0);
    }

    fprintf(stderr, "Galaxy radius read in.\n");

    // Import 4000A break (Dn4000).

    ff = "/Users/mehmet/Dropbox/nyucosmo/vagc_files_groups/dn4k.dr72bright34.dat";

    fp = fopen(ff,"r");
	if(!(fp=fopen(ff,"r"))){
      fprintf(stderr,"ERROR opening [%s]\n",ff);
      exit(0);
    } 

    j = 1;
    for(i = 1; i <= ngal; ++i){
    	if(mag_r[i] <= MAGNITUDE && redshift[i] <= MAXREDSHIFT && redshift[i] >= MINREDSHIFT && m_stellar[i] >= MSTARLIM){
    		fscanf(fp, "%f %*f %*d %*d %*f %*f %*f %*f %*f", &Dn4000[j]);
    		++j;
    	}
    }
	fclose(fp);

    if(j-1 != nsample){
    	printf("Error! Dn4000 array size not the same as sample size: %d %d\n", j-1, nsample);
    	exit(0);
    }

    fprintf(stderr, "Dn4000 read in.\n");

    // Import specific star formation rate (sSFR).

    ff = "/Users/mehmet/Dropbox/nyucosmo/vagc_files_groups/sfr.dr72bright34.dat";

    fp = fopen(ff,"r");
	if(!(fp=fopen(ff,"r"))){
    	fprintf(stderr,"ERROR opening [%s]\n",ff);
    	exit(0);
    }
    j = 1;
    for(i = 1; i <= ngal; ++i){
    	if(mag_r[i] <= MAGNITUDE && redshift[i] <= MAXREDSHIFT && redshift[i] >= MINREDSHIFT && m_stellar[i] >= MSTARLIM){
    		fscanf(fp, "%f %*f %*d %*d %*f %*f %*f %*f %*f", &sSFR[j]);
    		++j;
    	}
    }
    fclose(fp);

	if(j-1 != nsample){
    	printf("Error! sSFR array size not the same as sample size: %d %d\n", j-1, nsample);
    	exit(0);
    }

    fprintf(stderr, "sSFR read in.\n");

    // Import velocity dispersion and signal-to-noise (velDisp, sNr).

    ff = "/Users/mehmet/Dropbox/nyucosmo/vagc_files_groups/vdisp.dr72bright34.dat";

    fp = fopen(ff,"r");
	if(!(fp=fopen(ff,"r"))){
    	fprintf(stderr,"ERROR opening [%s]\n",ff);
    	exit(0);
    } 

    j =1;
    for(i = 1; i <= ngal; ++i){
    	if(mag_r[i] <= MAGNITUDE && redshift[i] <= MAXREDSHIFT && redshift[i] >= MINREDSHIFT && m_stellar[i] >= MSTARLIM){
    		fscanf(fp, "%f %*f %f", &velDisp[j], &sNr[j]);
    		++j;
    	}
    }
    fclose(fp);

    if(j-1 != nsample){
    	printf("Error! velDisp and/or sNr array size not the same as sample size: %d %d\n", j-1, nsample);
    	exit(0);
    }

    fprintf(stderr, "velDisp and sNr read in.\n");

    // Import Sersic index (sersic).

    ff = "/Users/mehmet/Dropbox/nyucosmo/vagc_files_groups/sersic_n.dr72bright34.dat";

    fp = fopen(ff,"r");
	if(!(fp=fopen(ff,"r"))){
    	fprintf(stderr,"ERROR opening [%s]\n",ff);
    	exit(0);
    }

    j = 1;
    for(i = 1; i <= ngal; ++i){
    	if(mag_r[i] <= MAGNITUDE && redshift[i] <= MAXREDSHIFT && redshift[i] >= MINREDSHIFT && m_stellar[i] >= MSTARLIM){
    		fscanf(fp, "%*d %f %*f %*f", &sersic[j]);
    		++j;
    	}
    }
    fclose(fp);

    if(j-1 != nsample){
    	printf("Error! sersic array size not the same as sample size: %d %d\n", j-1, nsample);
    	exit(0);
    }

    fprintf(stderr, "sersic read in.\n");

    // Import fibre collision (collision).

    ff = "/Users/mehmet/Dropbox/nyucosmo/vagc_files_groups/collided.dr72bright34.dat";

    fp = fopen(ff,"r");
	if(!(fp=fopen(ff,"r"))){
    	fprintf(stderr,"ERROR opening [%s]\n",ff);
    	exit(0);
    }

    j = 1;
    for(i = 1; i <= ngal; ++i){
    	if(mag_r[i] <= MAGNITUDE && redshift[i] <= MAXREDSHIFT && redshift[i] >= MINREDSHIFT && m_stellar[i] >= MSTARLIM){
    		fscanf(fp, "%d %*d", &collision[j]);
    		++j;
    	}
    }
    fclose(fp);

    if(j-1 != nsample){
    	printf("Error! collision array size not the same as sample size: %d %d\n", j-1, nsample);
    	exit(0);
    }

    fprintf(stderr, "collision read in.\n");

    // Import Petrosian radii (petroRad).

    ff = "/Users/mehmet/Dropbox/nyucosmo/vagc_files_groups/petrorad.dr72bright34.dat";

    fp = fopen(ff,"r");
	if(!(fp=fopen(ff,"r"))){
    	fprintf(stderr,"ERROR opening [%s]\n",ff);
    	exit(0);
    }

    j = 1;
    for(i = 1; i <= ngal; ++i){
    	if(mag_r[i] <= MAGNITUDE && redshift[i] <= MAXREDSHIFT && redshift[i] >= MINREDSHIFT && m_stellar[i] >= MSTARLIM){
    		fscanf(fp, "%f %f %*f %*f %*f", &x1, &x2);
    		petroRad[j] = x2/x1;
    		++j;
    	}
    }
    fclose(fp);

    if(j-1 != nsample){
    	printf("Error! petroRad array size not the same as sample size: %d %d\n", j-1, nsample);
    	exit(0);
    }

    fprintf(stderr, "petroRad read in.\n");

    // Before proceeding, go back through and truncate existing arrays to only contain galaxies from the volume limited sample.

  	float *temp_mag_r,*temp_mag_g,*temp_ra,*temp_dec,*temp_redshift,*temp_m_stellar;
  	int *temp_indx;

  	temp_mag_r = vector(1, nsample);
	temp_mag_g = vector(1, nsample);
	temp_ra = vector(1, nsample);
	temp_dec = vector(1, nsample);
	temp_redshift = vector(1, nsample);
	temp_m_stellar = vector(1, nsample);
	temp_indx = ivector(1, nsample);

  	j = 1;
  	for(i = 1; i <= ngal; ++i){
  		if(mag_r[i] <= MAGNITUDE && redshift[i] <= MAXREDSHIFT && redshift[i] >= MINREDSHIFT && m_stellar[i] >= MSTARLIM){
  			temp_mag_r[j] = mag_r[i];
  			temp_mag_g[j] = mag_g[i];
  			temp_ra[j] = ra[i];
  			temp_dec[j] = dec[i];
  			temp_redshift[j] = redshift[i];
  			temp_m_stellar[j] = m_stellar[i];
  			temp_indx[j] = j;
  			++j;
  		}
  	}

  	free(ra);
  	free(dec);
  	free(redshift);
  	free(mag_r);
  	free(mag_g);
  	free(m_stellar);
  	free(indx);

    ra = vector(1,nsample);
  	dec = vector(1,nsample);
  	redshift = vector(1,nsample);
  	mag_r = vector(1, nsample);
  	mag_g = vector(1, nsample);
  	m_stellar = vector(1, nsample);
  	indx = ivector(1, nsample);
  	indx = temp_indx;
  	ra = temp_ra;
  	dec = temp_dec;
  	redshift = temp_redshift;
  	for(i = 1; i <= nsample; ++i)
  		redshift[i] *= speedOfLight;
  	mag_r = temp_mag_r;
  	mag_g = temp_mag_g;
  	m_stellar = temp_m_stellar;

  	// Whew, all done. Now clear variables and proceed. This whole process should result in a series of arrays of length nsample that only contain galaxies in a volume limited sample.

    // Finished importing data!

    fprintf(stderr, "** Finished reading in data! Sample defined. **\n\n");
    printf("%d galaxies with mag_r > %2.2f, %.5f < z < %.4f, and mstar > %3.2e.\n",nsample,MAGNITUDE,MINREDSHIFT,MAXREDSHIFT, MSTARLIM);
    printf("Vol %e\n\n",volume);

    // Correct sSFR, Dn4000, and Hdelta for fibre collisions. Check fiber_corrected_galaxy_property.c for details.

    // This is in Jeremy's code but I'm not actually sure any of these changes are propagated into the rest of the analysis. Comment out for now.

	// for(i = 1; i <= ngal; ++i){
	// 	x1 = Dn4000[i];
	// 	x2 = sSFR[i];
	// 	x3 = Hdelta[i];
	// 	if(fabs(Dn4000[i] - 1.25) < 1.0E-4)
 //    		x1 = fiber_corrected_galaxy_property(i,mag_r,mag_g,Dn4000,ngal,collision);
 //    	if(fabs(Dn4000[i]-1.80)<1.0E-4)
	// 		x1 = fiber_corrected_galaxy_property(i,mag_r,mag_g,Dn4000,ngal,collision);
	// 	if(sSFR[i]<-98)
	// 		x3 = fiber_corrected_galaxy_property(i,mag_r,mag_g,sSFR,ngal,collision);
	// 	if(fabs(Hdelta[i])<1.0E-5)
	// 		x2 = fiber_corrected_galaxy_property(i,mag_r,mag_g,Hdelta,ngal,collision);
		
	// 	printf("GALDAT_CORR %d %f %f %f %f %f %f %e %f %f %f %f %f %f %d %f\n", i, mag_r[i], mag_g[i], redshift[i], x1, x2, x3, m_stellar[i], ra[i], dec[i], velDisp[i], sNr[i], sersic[i], petroRad[i],ka[i], Rexp[i]);
	// 	fflush(stdout);
 //    }

    // Sort everything by descending galaxy mass. In other words, run sort2 for every parameter that has been read in so far.

    for(i = 1; i <= nsample; ++i) {
    	m_stellar[i] = -m_stellar[i];
    	tempArray[i] = m_stellar[i];
    }
    // Store m_stellar in this temporary array. Before each sort, reassign m_stellar to this array; otherwise sorting won't happen as m_stellar will have already been sorted! Don't you just love C?

    sort2(nsample, tempArray, indx);
    for(i = 1; i <= nsample; ++i) tempArray[i] = m_stellar[i];
    sort3(nsample, tempArray, ra);
    for(i = 1; i <= nsample; ++i) tempArray[i] = m_stellar[i];
    sort3(nsample, tempArray, dec);
    for(i = 1; i <= nsample; ++i) tempArray[i] = m_stellar[i];
    sort3(nsample, tempArray, redshift);
    for(i = 1; i <= nsample; ++i) tempArray[i] = m_stellar[i];
    sort3(nsample, tempArray, mag_r);
    for(i = 1; i <= nsample; ++i) tempArray[i] = m_stellar[i];
    sort3(nsample, tempArray, mag_g);
    for(i = 1; i <= nsample; ++i) tempArray[i] = m_stellar[i];
    sort3(nsample, tempArray, v_max);
    for(i = 1; i <= nsample; ++i) tempArray[i] = m_stellar[i];
    sort3(nsample, tempArray, luminosity);
    for(i = 1; i <= nsample; ++i) tempArray[i] = m_stellar[i];
    sort3(nsample, tempArray, Hdelta);
    for(i = 1; i <= nsample; ++i) tempArray[i] = m_stellar[i];
    sort2(nsample, tempArray, ka);
    for(i = 1; i <= nsample; ++i) tempArray[i] = m_stellar[i];
    sort3(nsample, tempArray, Rexp);
    for(i = 1; i <= nsample; ++i) tempArray[i] = m_stellar[i];
    sort3(nsample, tempArray, Dn4000);
    for(i = 1; i <= nsample; ++i) tempArray[i] = m_stellar[i];
    sort3(nsample, tempArray, sSFR);
    for(i = 1; i <= nsample; ++i) tempArray[i] = m_stellar[i];
    sort3(nsample, tempArray, sersic);
    for(i = 1; i <= nsample; ++i) tempArray[i] = m_stellar[i];
    sort2(nsample, tempArray, collision);
    for(i = 1; i <= nsample; ++i) tempArray[i] = m_stellar[i];
    sort3(nsample, tempArray, velDisp);
    for(i = 1; i <= nsample; ++i) tempArray[i] = m_stellar[i];
    sort3(nsample, tempArray, sNr);
    for(i = 1; i <= nsample; ++i) tempArray[i] = m_stellar[i];
    sort3(nsample, m_stellar, petroRad);

    // m_stellar is now in descending order, with corresponding indx values. Make m_stellar positive again.

    for(i = 1; i <= nsample; ++i){
    	m_stellar[i] = -m_stellar[i];
    	luminosity[i] = m_stellar[i];
    }

    fprintf(stderr,"** Sorting complete! **\n\n");

    // ** UNCOMMENT BELOW TO OUTPUT A FILE WITH ALL THE DATA THAT HAS BEEN READ IN **

   //  fp = fopen("/Users/mehmet/Desktop/v2in.dat","w");
  	// fprintf(fp,"ra,dec,redshift,mag_r,mag_g,v_max,m_stellar,Hdelta,ka,Rexp,Dn4000,sSFR,sersic,collision,velDisp,petroRad\n");
  	// for(i = 1; i <= nsample; ++i){
   //  	fprintf(fp,"%f,%f,%f,%f,%f,%f,%f,%f,%d,%f,%f,%f,%f,%d,%f,%f\n", ra[i],dec[i],redshift[i], mag_r[i], mag_g[i], v_max[i],m_stellar[i],Hdelta[i],ka[i],Rexp[i],Dn4000[i],sSFR[i],sersic[i],collision[i],velDisp[i],petroRad[i]);
  	// }
  	// fclose(fp);

    // ** UNCOMMENT ABOVE TO OUTPUT A FILE WITH ALL THE DATA THAT HAS BEEN READ IN **

	// Cycle through each galaxy and assign a halo mass to it. This part of the code calls up the density2halo function, which does SHAM to assign a halo mass to each galaxy.

	// This is the initial first pass assignment of masses to galaxies

    j = 0;
    k = 0;

    for(j = 1; j <= nsample; ++j){
    	i = indx[j];
    	ndens_gal += 1/v_max[i];
    	mass[i] = density2halo(ndens_gal);
    	rad[i] = pow((3*mass[i]) / (4.0 * pi * dHalo *rhoCrit * omegaM), third);
    	angRad[i] = rad[i] / distance_redshift(redshift[i]/speedOfLight);
    	sigma[i] = sqrt((bigG*mass[i])/(2.0*rad[i]));
    }
    
	fprintf(stderr, "** Initial assignment of halo masses to galaxies complete. **\n\n");

	//fprintf(stderr, "Number density = %3.3e\n\n",ndens_gal);

	// Now go through and find associated galaxies.

	fprintf(stderr, "** Identifying satellites...\n\n");

	// Prepare variables/arrays.
	// igrp is simply the group name. 
	// group_member keeps track of which galaxy belongs in which group name.
	// group_center keeps track of which galaxy is at the centre of each group.

	k = 0;
	igrp = 0;
	nsat_tot = 0;
	for(i = 1; i <= nsample; ++i){
		
		prob_total[i] = 0.;
		group_member[i] = 0;
	}
	
	for(i = 1; i <= nsample; ++i){
		j = indx[i];
		igrp++;
		group_luminosity[igrp] = m_stellar[i];
		group_center[igrp] = j;
		group_member[j] = igrp; 

		// if(prob_total[j] > 0)
			// printf("%d %d %f\n",i,j,prob_total[j]);

		group_luminosity[igrp] += find_satellites(j, ra, dec, redshift, mag_r, angRad[j], sigma[j], group_member, indx, nsample, rad[j], mass[j], igrp, luminosity,&nsat_indi[j],i,prob_total);
		
		if(nsat_indi[j] < 1)
			k++;
		nsat_tot += nsat_indi[j] * (1-prob_total[j]);
		
	}

	fprintf(stderr, "** First pass satellite identification complete: \n");
	ngrp = igrp;
	
	printf("%d groups, %d n = 1 groups, fsat = %.2f **\n\n", ngrp, k, nsat_tot*1./nsample);

	// We now have the first iteration of groups. As before, rank these by their group mass.

	for(i = 1; i <= ngrp; ++i)
		group_luminosity[i] *= -1;
	sort2(ngrp, group_luminosity, group_center);
	for(i = 1; i <= ngrp; ++i)
		group_luminosity[i] *= -1;

	// Now repeat the abundance matching, but this time use total group mass instead.

	j = 0;
	xngrp = 0;
	for(j = 1; j <= ngrp; ++j){
		i = group_center[j];
		mass[i] = density2host_halo(j/volume);
		rad[i] = pow(3*mass[i]/(4.*pi*dHalo*rhoCrit*omegaM),1.0/3.0);
		angRad[i] = rad[i]/distance_redshift(redshift[i]/speedOfLight);
		sigma[i] = sqrt((bigG*mass[i])/(2.0*rad[i])*(1+redshift[i]/speedOfLight));
	}

	fprintf(stderr, "** SHAMmed groups! Now iterating to convergence... ** \n\n");

	// This is the main iteration loop.

	for(niter = 1; niter <= niter_max; ++niter){

		// Begin by resetting group membership.

		printf("Iteration %d of %d...\n", niter, niter_max);

		j = 0;
		for(i = 1; i <= nsample; ++i){
			if(prob_total[i] >= 0)
				prob_total[i] = 0.;
		}
		for(i = 1; i <= nsample; ++i){
			group_member[i] = 0;
			nsat_indi[i] = 0;
		}

		// Transfer group centers to a temporary list.

		for(i = 1; i <= ngrp; ++i){
			temp_group[i] = group_center[i];
		}
		ngrp_temp = ngrp;
		ngrp = 0;
		k = 0;
		nsat_tot = 0;

		// Reset complete!

		// Perform satellite identification for current iteration of the groups.

		for(j = 1; j <= ngrp_temp; ++j){

			i = temp_group[j];

			// Is this galaxy a group member? If so, skip!

			if(group_member[i] > 0)
				continue;
			ngrp++;
			group_member[i] = ngrp;
			group_center[ngrp] = i;
			nsat_indi[i] = 0;

			group_luminosity[ngrp] += find_satellites(i, ra, dec, redshift, mag_r, angRad[i], sigma[i], group_member, indx, nsample, rad[i], mass[i], ngrp, luminosity, &nsat_indi[i], j, prob_total);

			if(nsat_indi[i] == 0)
				k++;

			nsat_tot += nsat_indi[i];

		}

		// Some galaxies will now have been newly 'exposed.'

		for(i = 1; i <= nsample; ++i){
			
			if(group_member[i])
				continue;

			ngrp++;
			group_luminosity[ngrp] = m_stellar[i];
			group_member[i] = ngrp;
			group_center[ngrp] = i;

			nsat_indi[i] = 0;
			group_luminosity[ngrp] += find_satellites(i, ra, dec, redshift, mag_r, angRad[i], sigma[i], group_member, indx, nsample, rad[i], mass[i], ngrp, luminosity, &nsat_indi[i], j, prob_total);

			if(nsat_indi[i] == 0)
				k++;
			nsat_tot += nsat_indi[i];
		}

		printf("%d groups, %d n = 1 groups, fsat = %.2f\n", ngrp, k, nsat_tot*1./nsample);

	}

}

// ** End of main program. **

// ** Local functions. **

float func_dr1(float z) {
	return pow(omegaM*(1+z)*(1+z)*(1+z)+(1-omegaM),-0.5);
}

float distance_redshift(float z) {
	float x;
	if(z<=0)
		return 0;
	x = cOnH0 * qromo(func_dr1, 0.0, z, midpnt);
	return x;
}

float angular_separation(float a1, float d1, float a2, float d2)
{
 // float cd1,cd2,sd1,sd2,ca1a2,sa1a2;

  return atan((sqrt(cos(d2)*cos(d2)*sin(a2-a1)*sin(a2-a1) + 
		    pow(cos(d1)*sin(d2) - sin(d1)*cos(d2)*cos(a2-a1),2.0)))/
	      (sin(d1)*sin(d2) + cos(d1)*cos(d2)*cos(a2-a1)));
}

float radial_probability(float mass, float dr, float rad, float ang_rad)
{
  float c, x, rs, delta, f;

  dr = dr*rad/ang_rad;

  c = 10.0*pow(mass/1.0E+14,-0.11);
  rs = rad/c;
  x = dr/rs;

  if(x<1)
    f = 1/(x*x-1)*(1-log((1+sqrt(1-x*x))/x)/(sqrt(1-x*x)));
  if(x==1)
    f = 1.0/3.0;
  if(x>1)
    f = 1/(x*x-1)*(1-atan(sqrt(x*x-1))/sqrt(x*x-1));

  delta = dHalo/3.0*c*c*c/(log(1+c)-c/(1+c));

  return 1.0/cOnH0*2*rs*delta*f;
}

// find_satellites identifies satellite galaxies near the central galaxy (whose coordinates are supplied in the argument to the function). 

float find_satellites(int i, float *ra, float *dec, float *redshift, float *mag_r, float theta_max, float sigma, int *group_member, int *indx, int ngal, float radius, float mass, int igrp, float *luminosity, float *nsat_cur, int i1, float *prob_total) {
	int j, k;
	float dx, dy, dz, theta, prob_ang, vol_corr, prob_rad, grp_lum, p0;
	
	// Set up.

	*nsat_cur = 0;
	grp_lum = 0;

	for(k = 1; k <= ngal; ++k){
		j = indx[k];

		// The idea is to loop through every other galaxy in the sample and calculate the probability that it's in the same halo. This normally takes ages, we include some checks.

		// Skip if target galaxy is the same as the central (obviously).
		if(j == i)
		 	continue;

		// Skip if already assigned to a central.
		if(group_member[j]){
			continue;
		}
		
		// This chunk of code truncates the full sample of galaxies to a narrow 'window' around the target central galaxy in order to reduce computation time. This will be replaced with a k-d tree search.

		// ** INJECT K-D TREE BELOW ** //
		
		dx = fabs(ra[i] - ra[j]);
		if(dx > 4*theta_max)
			continue;
		dy = fabs(dec[i] - dec[j]);
		if(dy > 4*theta_max)
			continue;
		dz = fabs(redshift[i] - redshift[j]);
		if(dz > 6*sigma)
			continue;
				
		// ** INJECT K-D TREE ABOVE ** //

		theta = angular_separation(ra[i],dec[i],ra[j],dec[j]);

		if(theta > theta_max)
			continue;

		// Now determine the probability of being a satellite (both projected onto the sky, and along the line of sight).

		prob_ang = radial_probability(mass, theta, radius, theta_max);
		prob_rad = (exp(-dz*dz/(2*sigma*sigma)) * speedOfLight) / (rt2Pi*sigma);
		
		p0 = (1 - 1/(1 + prob_ang * prob_rad / 10));
		
		if(p0 < 0){
			printf("ZERO %e\n",p0);
			p0 = 0;
		}

		if(p0 > prob_total[j])
			prob_total[j] = p0;
		if(prob_total[j] > 1)
			prob_total[j] = 1;
		if(prob_ang*prob_rad < 10)
			continue;
		
		// At this point the galaxy is a satellite. Assign the group ID number to it, and add its mass to the total group mass. Increase satellite counter by 1.

		group_member[j] = igrp;
		grp_lum += luminosity[k];
		(*nsat_cur) += 1;

	}
	
	// Vmax correction.

	dz = speedOfLight* fabs(redshift[i] - MINREDSHIFT);
	vol_corr = 1-(0.5*erfc(dz/(root2*sigma)));
	*nsat_cur /= vol_corr;
	grp_lum /= vol_corr;
	
	dz = speedOfLight* fabs(redshift[i] - MAXREDSHIFT);
	vol_corr = 1-(0.5*erfc(dz/(root2*sigma)));
	*nsat_cur /= vol_corr;
	grp_lum /= vol_corr;

	return(grp_lum);

}