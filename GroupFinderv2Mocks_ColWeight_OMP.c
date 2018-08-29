// As GroupFinderv2, but designed to run on mock catalogue. Modified to include a weighting factor to upscale or downscale stellar masses of red centrals.

// For self: compile command --
// gcc -o gfv2mockColOMP GroupFinderv2Mocks_ColWeight_OMP.c *.o -L/home/users/ma5046/libC_main -lC_main -lm -fopenmp
//
// Initialization //

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <sys/time.h>
#include "nrutil.h"
#include "kdtree.h"
#include "omp.h"

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
#define czMax 25502.0
#define czMin 6001.0
#define czBuf 0
#define CHUNKSIZE 10000
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
float find_satellites(int i, float *ra, float *dec, float *redshift, int *color_flag, float theta_max, float x1, int *group_member, int *indx, int ngal, float radius, float mass, int igrp, float *m_stellar, float *nsat_cur, float *prob_total, void *kd, float cenDist, float range, float sw1, float sw2);
float radial_probability(float mass, float dr, float rad, float ang_rad);
int iter_central_galaxy(int galID, float *ra, float *dec, float *redshift ,int *group_member, int ngal, int igrp, float *m_stellar, float *massDist);
float segvol(float lorad, float hirad, float lodec, float hidec, float lora, float hira);
float segarea(float lorad, float hirad, float lodec, float hidec, float lora, float hira);

//static double dist_sq( double *a1, double *a2, int dims );

unsigned int get_msec(void)
{
  static struct timeval timeval, first_timeval;

  gettimeofday(&timeval, 0);

  if(first_timeval.tv_sec == 0) {
    first_timeval = timeval;
    return 0;
  }
  return (timeval.tv_sec - first_timeval.tv_sec) * 1000 + (timeval.tv_usec - first_timeval.tv_usec) / 1000;
}

// Global variables go here

float GALAXY_DENSITY, MAGNITUDE, MAXREDSHIFT, MINREDSHIFT;

// End initialization //

// Here we go!

int main(int argc, char **argv){

  int i, j, k, igrp, ngal, nsample, count, imax, *GalID;
  int *indx, *HaloID, *central_flag, *color_flag;
  int *group_member, *group_index, *group_center, *temp_group;
  int nsat_tot, ngrp, niter, niter_max, ngrp_temp;
  int chunk = CHUNKSIZE;

  float *ra, *dec, *redshift, *m_stellar, *m_halo;
  float *mass, *rad, *angRad, *sigma, *prob_total, *nsat, maxMass, *cenDist, *range, sw1, sw2;
  float *group_mass, *distToBig;
  float x1, x2, *tempArray;
  float volume;
  float red_weight;
  float x, y, z, radius, theta;
  unsigned int msec, start;

  void *kd;

  double ndens_gal = 0;

  char string[1000];
  char *ff;
  FILE *fp;

  clock_t startinit, endinit, start1, end1, start2, end2, start3, end3, startout, endout, startall, endall;
  float *zone1, *zone2, *zone3, zone1t, zone2t, zone3t, zoneinit, zoneout;

  startall = get_msec();
  startinit = get_msec();
  count = 0;
  niter_max = 10;

  zone1 = vector(1, niter_max);
  zone2 = vector(1, niter_max);
  zone3 = vector(1, niter_max);

  MAXREDSHIFT = czMax / speedOfLight;
  MINREDSHIFT = czMin / speedOfLight;
  MAGNITUDE = -19;
  //MSTARLIM = pow(10.0,15);
  x1 = x2 = 0;

  volume = segvol(distance_redshift(MINREDSHIFT),distance_redshift(MAXREDSHIFT), 0, 90, 0, 90);
  
  printf("Magitude: %f\nMin redshift: %f\nMax redshift: %f\nVolume: %e\n",MAGNITUDE,MINREDSHIFT,MAXREDSHIFT,volume);

  // Start by reading in data. At first only read in bare minimum to establish a volume limited sample; then read in everything else.

  printf("\n** Reading in data and defining sample. **\n");

  // Import mock catalogue (RA, Dec, z). 
  // Measure length of this catalogue; this is the total number of galaxies. Apply this length to the arrays containing redshift, magnitude, and mass.
  //
  ff = "/home/users/ma5046/misc_work/mocks/bolshoiColor_input.csv";
  
  fp = fopen(ff,"r");
  if(!(fp=fopen(ff,"r"))){
      printf("ERROR opening [%s]\n",ff);
      exit(0);
    }

  // Skip first line of this file. It will have a header line, as it is a CSV.
  fgets(string,1000,fp);

  // This measures the length of the file.

  while(fgets(string,1000,fp)){
      count++;
  }
  rewind(fp);
  // Assign count to ngal; this is the number of galaxies.

  ngal = count;
  printf("Red weight is %f\n",atof(argv[1]));
  // Each variable that holds a different galaxy property can now be made into an array going from 1 to ngal.

  ra = vector(1,ngal);
  dec = vector(1,ngal);
  redshift = vector(1,ngal);
  m_stellar = vector(1, ngal);
  m_halo = vector(1, ngal);
  central_flag = ivector(1, ngal);
  color_flag = ivector(1, ngal);
  indx = ivector(1, ngal);
  HaloID = ivector(1, ngal);
  GalID = ivector(1, ngal);

  for(i = 1; i <= ngal; ++i){

    // Skip first line of mocks_HaloID.csv as this contains header info.

    if(i == 1)
      fgets(string,1000,fp);

    fscanf(fp,"%f,%f,%f,%f,%f,%d,%d,%d",&ra[i],&dec[i],&redshift[i], &m_stellar[i], &m_halo[i], &HaloID[i], &central_flag[i], &color_flag[i]);
    GalID[i] = i;
    ra[i] *= pi/180;
    dec[i] *= pi/180;
    indx[i] = i;
    m_stellar[i] = pow(10.0,m_stellar[i]);
    fgets(string,1000,fp);
  }

  fclose(fp);
  printf("Mock data read in.\n");

  // Now have enough data to make a volume limited sample.
  // Start by identifying how many galaxies are in this sample. Assign this length to all remaining data arrays.

  count = 0;
  for(i = 1; i <= ngal; ++i){
    if(redshift[i] <= MAXREDSHIFT && redshift[i] >= MINREDSHIFT){
      count++;
    }
  }

  nsample = count;
  printf("%d %d\n",ngal,nsample);

  mass = vector(1, nsample);
  rad = vector(1, nsample);
  angRad = vector(1, nsample);
  sigma = vector(1, nsample);
  prob_total = vector(1, nsample);

  tempArray = vector(1, nsample);
  temp_group = ivector(1, nsample);

  group_mass = vector(1, nsample);
  nsat = vector(1, nsample);
  distToBig = vector(1, nsample);
  group_member = ivector(1, nsample);
  group_index = ivector(1, nsample);
  group_center = ivector(1, nsample);

  // Before proceeding, go back through and truncate existing arrays to only contain galaxies from the volume limited sample.

  float *temp_ra,*temp_dec,*temp_redshift,*temp_m_stellar, *temp_m_halo;
  int *temp_indx, *temp_HaloID, *temp_central_flag, *temp_origGalID, *temp_color_flag;

  temp_ra = vector(1, nsample);
  temp_dec = vector(1, nsample);
  temp_redshift = vector(1, nsample);
  temp_m_stellar = vector(1, nsample);
  temp_m_halo = vector(1, nsample);
  temp_indx = ivector(1, nsample);
  temp_HaloID = ivector(1, nsample);
  temp_origGalID = ivector(1, nsample);
  temp_central_flag = ivector(1, nsample);
  temp_color_flag = ivector(1, nsample);

    j = 1;
    for(i = 1; i <= ngal; ++i){
      if(redshift[i] <= MAXREDSHIFT && redshift[i] >= MINREDSHIFT){
        temp_ra[j] = ra[i];
        temp_dec[j] = dec[i];
        temp_redshift[j] = redshift[i];
        temp_m_stellar[j] = m_stellar[i];
        temp_m_halo[j] = m_halo[i];
        temp_HaloID[j] = HaloID[i];
        temp_central_flag[j] = central_flag[i];
        temp_color_flag[j] = color_flag[i];
        temp_origGalID[j] = GalID[i];
        temp_indx[j] = j;
        ++j;
      }
    }

    free(ra);
    free(dec);
    free(redshift);
    free(m_stellar);
    free(m_halo);
    free(HaloID);
    free(central_flag);
    free(color_flag);
    free(indx);
    free(GalID);

    ra = vector(1,nsample);
    dec = vector(1,nsample);
    redshift = vector(1,nsample);
    m_stellar = vector(1, nsample);
    m_halo = vector(1, nsample);
    HaloID = ivector(1, nsample);
    central_flag = ivector(1, nsample);
    color_flag = ivector(1, nsample);
    indx = ivector(1, nsample);
    GalID = ivector(1,nsample);

    indx = temp_indx;
    ra = temp_ra;
    dec = temp_dec;
    redshift = temp_redshift;
    for(i = 1; i <= nsample; ++i)
      redshift[i] *= speedOfLight;
    m_stellar = temp_m_stellar;
    m_halo = temp_m_halo;
    HaloID = temp_HaloID;
    GalID = temp_origGalID;
    central_flag = temp_central_flag;
    color_flag = temp_color_flag;

    // Whew, all done. Now clear variables and proceed. This whole process should result in a series of arrays of length nsample that only contain galaxies in a volume limited sample.

    // Finished importing data!

    printf("** Finished reading in data! Sample defined. **\n\n");
    printf("%d galaxies with %.5f < z < %.5f.\n",nsample,MINREDSHIFT,MAXREDSHIFT);
    printf("Vol %e\n\n",volume);

    // Correct sSFR, Dn4000, and Hdelta for fibre collisions. Check fiber_corrected_galaxy_property.c for details.

    // This is in Jeremy's code but I'm not actually sure any of these changes are propagated into the rest of the analysis. Comment out for now.

  // for(i = 1; i <= ngal; ++i){
  //  x1 = Dn4000[i];
  //  x2 = sSFR[i];
  //  x3 = Hdelta[i];
  //  if(fabs(Dn4000[i] - 1.25) < 1.0E-4)
 //       x1 = fiber_corrected_galaxy_property(i,mag_r,mag_g,Dn4000,ngal,collision);
 //     if(fabs(Dn4000[i]-1.80)<1.0E-4)
  //    x1 = fiber_corrected_galaxy_property(i,mag_r,mag_g,Dn4000,ngal,collision);
  //  if(sSFR[i]<-98)
  //    x3 = fiber_corrected_galaxy_property(i,mag_r,mag_g,sSFR,ngal,collision);
  //  if(fabs(Hdelta[i])<1.0E-5)
  //    x2 = fiber_corrected_galaxy_property(i,mag_r,mag_g,Hdelta,ngal,collision);
    
  //  printf("GALDAT_CORR %d %f %f %f %f %f %f %e %f %f %f %f %f %f %d %f\n", i, mag_r[i], mag_g[i], redshift[i], x1, x2, x3, m_stellar[i], ra[i], dec[i], velDisp[i], sNr[i], sersic[i], petroRad[i],ka[i], Rexp[i]);
  //  fflush(stdout);
 //    }

    // Sort everything by descending galaxy mass. In other words, run sort2 for every parameter that has been read in so far.

    for(i = 1; i <= nsample; ++i) {
      m_stellar[i] = -m_stellar[i];
      tempArray[i] = m_stellar[i];
    }
    // Store m_stellar in this temporary array. Before each sort, reassign m_stellar to this array; otherwise sorting won't happen as m_stellar will have already been sorted! Don't you just love C? Make sure the last one isn't temparray though, or m_stellar won't get sorted.


    sort2(nsample, tempArray, indx);
    for(i = 1; i <= nsample; ++i) tempArray[i] = m_stellar[i];
    sort3(nsample, tempArray, ra);
    for(i = 1; i <= nsample; ++i) tempArray[i] = m_stellar[i];
    sort3(nsample, tempArray, dec);
    for(i = 1; i <= nsample; ++i) tempArray[i] = m_stellar[i];
    sort3(nsample, tempArray, redshift);
    for(i = 1; i <= nsample; ++i) tempArray[i] = m_stellar[i];
    sort3(nsample, tempArray, m_halo);
    for(i = 1; i <= nsample; ++i) tempArray[i] = m_stellar[i];
    sort2(nsample, tempArray, HaloID);
    for(i = 1; i <= nsample; ++i) tempArray[i] = m_stellar[i];
    sort2(nsample, tempArray, central_flag);
    for(i = 1; i <= nsample; ++i) tempArray[i] = m_stellar[i];
    sort2(nsample, tempArray, color_flag);
    for(i = 1; i <= nsample; ++i) tempArray[i] = m_stellar[i];
    sort2(nsample, m_stellar, GalID);

    // m_stellar is now in descending order, with corresponding indx values. Make m_stellar positive again.

    for(i = 1; i <= nsample; ++i){
      m_stellar[i] = -m_stellar[i];
    }

    // While we're at it, compute distances to each galaxy too. This makes life easier later with OpenMP.
    
    cenDist = vector(1, nsample);
    range = vector(1, nsample);
    
    for(i = 1; i <= nsample; ++i){
      cenDist[i] = distance_redshift(redshift[i]/speedOfLight);
    }

    printf("** Sorting complete! **\n\n");

    // ** UNCOMMENT BELOW TO OUTPUT A FILE WITH ALL THE DATA THAT HAS BEEN READ IN **

   //  fp = fopen("/Users/mehmet/Desktop/v2in.dat","w");
    // fprintf(fp,"ra,dec,redshift,mag_r,mag_g,v_max,m_stellar,Hdelta,ka,Rexp,Dn4000,sSFR,sersic,collision,velDisp,petroRad\n");
    // for(i = 1; i <= nsample; ++i){
      // fprintf(fp,"%f,%f,%f,%f,%f,%f,%f,%f,%d,%f,%f,%f,%f,%d,%f,%f\n", ra[i],dec[i],redshift[i], mag_r[i], mag_g[i], v_max[i],m_stellar[i],Hdelta[i],ka[i],Rexp[i],Dn4000[i],sSFR[i],sersic[i],collision[i],velDisp[i],petroRad[i]);
    // }
    // fclose(fp);

    // ** UNCOMMENT ABOVE TO OUTPUT A FILE WITH ALL THE DATA THAT HAS BEEN READ IN **

  // Cycle through each galaxy and assign a halo mass to it. This part of the code calls up the density2halo function, which does SHAM to assign a halo mass to each galaxy.

  // This is the initial first pass assignment of masses to galaxies

    j = 0;
    k = 0;

    for(i = 1; i <= nsample; ++i){
      ndens_gal += 1/volume;
      mass[i] = density2halo(ndens_gal);
      rad[i] = pow((3*mass[i]) / (4.0 * pi * dHalo *rhoCrit * omegaM), third);
      angRad[i] = rad[i] / distance_redshift(redshift[i]/speedOfLight);
      sigma[i] = sqrt((bigG*mass[i])/(2.0*rad[i])*(1+redshift[i]/speedOfLight));
      range[i] = distance_redshift(4.0*sigma[i]/speedOfLight);
    }
    
  printf("** SHAMmed galaxies. **\n\n");

  //printf("Number density = %3.3e\n\n",ndens_gal);

  // Now go through and find associated galaxies.
  // Send satellite color arguments to these floats, to then feed to find_satellites.
  
  sw1 = atof(argv[4]);
  sw2 = atof(argv[3]);

  printf("** Identifying satellites...\n\n");

  // Create a k-d tree of all sample galaxies. This is used for nearest-neighbour searches later when identifying satellites around centrals.

    // k-d tree should have 3 dimensions (RA, Dec, z).

    kd = kd_create(3);

    // Insert points into tree. Each point consists of RA and Dec of a galaxy, projected onto a plane using the Hammer projection.

    // Names array is made to simply pass the row numbers of galaxies into the kd-tree for later retrieval.
    int *names;
    names = ivector(1, nsample);

    for(i = 1; i <= nsample; ++i){
      radius = distance_redshift(redshift[i]/speedOfLight);
      x = radius * cos(ra[i]) * cos(dec[i]);
      y = radius * sin(ra[i]) * cos(dec[i]); 
      z = radius * sin(dec[i]);
      names[i] = i;

      double pt[3] = {x, y, z};
      assert( kd_insert(kd, pt, (void*)&names[i]) == 0);
    }

  // Prepare variables/arrays.
  // igrp is simply the group name. 
  // group_member keeps track of which galaxy belongs in which group name.
  // group_center keeps track of which galaxy is at the centre of each group.
  // group_index is the same as indx but for groups.

  k = 0;
  igrp = 0;
  nsat_tot = 0;
  for(i = 1; i <= nsample; ++i){
    prob_total[i] = 0.;
    group_member[i] = 0;
  }

  //#pragma omp parallel num_threads (24) default(shared) private(i)
  //{
  //  #pragma omp for schedule (dynamic, chunk)
     
  //start = get_msec();
    for(i = 1; i <= nsample; ++i){
      
      if(group_member[i])
        continue;
      igrp++;
      
      // This is where we insert the weighting factor. If the central galaxy is red, apply the red weight to its stellar mass.

      if(color_flag[i] == 1) 
        red_weight = (atof(argv[2]) * log10(m_stellar[i])) + atof(argv[1]);
      else 
        red_weight = 1.0;
      group_mass[igrp] = m_stellar[i] * red_weight;
       
      group_center[igrp] = i;
      group_member[i] = igrp;

      group_mass[igrp] += find_satellites(i, ra, dec, redshift, color_flag, angRad[i], sigma[i], group_member, indx, nsample, rad[i], mass[i], igrp, m_stellar,&nsat[i],prob_total, kd, cenDist[i], range[i], sw1, sw2);
      if(nsat[i] < 1)
        k++;
      nsat_tot += nsat[i];
    }
  //}
  //msec = get_msec() - start;
  printf("** First pass satellite identification complete: \n");
  ngrp = igrp;
  
  //printf("%.3f sec\n", (float)msec / 1000.0);
  printf("%d groups, %d n = 1 groups, fsat = %.2f **\n\n", ngrp, k, nsat_tot*1./nsample);

  // We now have the first iteration of groups. As before, rank these by their group mass.

  for(i = 1; i <= ngrp; ++i)
    group_mass[i] *= -1;
  sort2(ngrp, group_mass, group_center);
  for(i = 1; i <= ngrp; ++i)
    group_mass[i] *= -1;

  // Now repeat the abundance matching, but this time use total group mass instead.

  for(i = 1; i <= ngrp; ++i){
    // skip if j i satellite, then go to ncount for density
    mass[i] = density2host_halo(i/volume);
    rad[i] = pow(3*mass[i]/(4.*pi*dHalo*rhoCrit*omegaM),1.0/3.0);
    angRad[i] = rad[i]/distance_redshift(redshift[i]/speedOfLight);
    sigma[i] = sqrt((bigG*mass[i])/(2.0*rad[i])*(1+redshift[i]/speedOfLight));
    range[i] = distance_redshift(4.0*sigma[i]/speedOfLight);
  }

  printf("** SHAMmed groups! Now iterating to convergence... ** \n\n");

  // This is the main iteration loop.
  endinit = get_msec() - startinit;
  zoneinit = (float)endinit / 1000.;

  for(niter = 1; niter <= niter_max; ++niter){
  
    // Begin by resetting group membership.

    start1 = get_msec();

    printf("Iteration %d of %d... \n", niter, niter_max);

    for(i = 1; i <= nsample; ++i){
      if(prob_total[i] >= 0)
        prob_total[i] = 0;
    }
    for(i = 1; i <= nsample; ++i){
      group_member[i] = 0;
      nsat[i] = 0;
    }

    // Transfer group centers to a temporary list.

    for(i = 1; i <= ngrp; ++i){
      temp_group[i] = group_center[i];
      distToBig[i] = 0.;
    }

    ngrp_temp = ngrp;
    igrp = ngrp = 0;
    k = 0;
    nsat_tot = 0;

    // Reset complete!

    end1 = get_msec() - start1;
    zone1[niter] = (float)end1 / 1000.;
    
    // Perform satellite identification for current iteration of the groups.

    start2 = get_msec();
    //#pragma omp parallel num_threads (24) default(shared) private(i)
    //  {
    //  #pragma omp for schedule (dynamic, chunk)    
      for(i = 1; i <= ngrp_temp; ++i){

        // i = temp_group[j];

        // Is this galaxy a group member? If so, skip!

        if(group_member[i]){
          continue;
        }

        igrp++;

        // Apply red weighting again.

        if(color_flag[i] == 1) 
          red_weight = (atof(argv[2]) * log10(m_stellar[i])) + atof(argv[1]);
        else 
          red_weight = 1.0;
        
        group_mass[igrp] = m_stellar[i] * red_weight;
        group_member[i] = igrp;
        group_center[igrp] = i;

        group_mass[igrp] += find_satellites(i, ra, dec, redshift, color_flag, angRad[i], sigma[i], group_member, indx, nsample, rad[i], mass[i], igrp, m_stellar, &nsat[i], prob_total, kd, cenDist[i], range[i], sw1, sw2);

        if(nsat[i] == 0)
          k++;

        nsat_tot += nsat[i];
      }
    
  
    // Some galaxies will now have been newly 'exposed.'

    count = 0;
    
      for(i = 1; i <= nsample; ++i){

        if(group_member[i]){
          continue;
        }

        if(prob_total[i] > 0.5){
          continue;
        }

        igrp++;

        // No red weighting this time, because we only want to do this for centrals.

        group_mass[igrp] = m_stellar[i];
        group_member[i] = igrp;
        group_center[igrp] = i;

        nsat[i] = 0;

        group_mass[igrp] += find_satellites(i, ra, dec, redshift, color_flag, angRad[i], sigma[i], group_member, indx, nsample, rad[i], mass[i], igrp, m_stellar, &nsat[i], prob_total, kd, cenDist[i], range[i], sw1, sw2);

        if(nsat[i] == 0)
          k++;
        nsat_tot += nsat[i];
      }
    //}
  
    ngrp = igrp;
    printf("%d groups, %d n = 1 groups, fsat = %.2f\n", ngrp, k, nsat_tot*1./nsample);
    for(i = 1; i <= ngrp; ++i){
      group_index[i] = i;
    }

    // Sort this iteration's groups by mass (descending).

    for(i = 1; i <= ngrp; ++i)
      group_mass[i] *= -1;
    sort2(ngrp, group_mass, group_index);
    for(i = 1; i <= ngrp; ++i)
      group_mass[i] *= -1;

    // Now perform abundance matching on this new set of groups. Must first determine new group centres based on updated group composition, however. Make sure not to run on satellites.
    end2 = get_msec() - start2;
    zone2[niter] = (float)end2 / 1000.;

    start3 = get_msec();
    #pragma omp parallel num_threads(24) default(shared) private(i,j) 
    {
    for(i = 1; i <= ngrp; ++i){

      igrp = group_index[i];
      k = group_center[igrp];

      // What's the new group center? Also compute distance between most massive galaxy and the group center.
      
      //if(nsat[k] > 3){
      //  j = iter_central_galaxy(k, ra, dec, redshift, group_member, nsample, igrp, m_stellar, &distToBig[i]);
    
      //  if(j != k){
      //    group_center[igrp] = j;
      //    k = j;
      //  }
      //}

      // New group centres have now been identified. Time to SHAM!

      mass[k] = density2host_halo(i/volume);
      rad[k] = pow(3*mass[k]/(4.*pi*dHalo*rhoCrit*omegaM),1.0/3.0);
      angRad[k] = rad[k]/distance_redshift(redshift[k]/speedOfLight);
      sigma[k] = sqrt((bigG*mass[k])/(2.0*rad[k])*(1+redshift[k]/speedOfLight));
      //range[k] = distance_redshift(4.0*sigma[k]/speedOfLight);
      //printf("BGH%d %f %f\n",niter,mass[k],group_mass[i]);
     
      // Identify most massive galaxy in group.
      
      maxMass = 0;


      #pragma omp for schedule (dynamic, chunk)
      for(j = 1; j <= nsample; ++j){
        if(group_member[j] == igrp){
          if(m_stellar[j] > maxMass){
            maxMass = m_stellar[j];
            imax = j;
          }
        }
      }
    }
    }
    end3 = get_msec() - start3;
    zone3[niter] = (float)end3 / 1000.;
    // If niter = niter_max, time to output!

    if(niter == niter_max){
      startout = get_msec();
      printf("** Iterations complete! Writing output to file... **\n\n");

      ff = "/home/users/ma5046/misc_work/output/bolshoiz01cGroupsOMP.csv";
      fp = fopen(ff,"w");
      fprintf(fp,"groupID,ra,dec,redshift,centralID,nsat,MSgroup,Mhalo,rad,sigma,angRad,Mcentral,massSep,HaloID,SimHaloMass\n");
      for(i = 1; i <= ngrp; ++i){
        j = group_index[i];
        k = group_center[j];
        fprintf(fp,"%d,%f,%f,%f,%d,%f,%f,%f,%f,%f,%f,%f,%f,%d,%f\n", j, ra[k], dec[k], redshift[k], k, nsat[k], group_mass[i], mass[k], rad[k], sigma[k], angRad[k],m_stellar[k],distToBig[k],HaloID[k],m_halo[k]);
      }
      fclose(fp);

      ff = "/home/users/ma5046/misc_work/output/bolshoiz01cGalsOMP.csv";
      fp = fopen(ff,"w");
      fprintf(fp,"galID,ra,dec,redshift,Mstellar,groupID,prob_total,centralID,Mhalo,Rhalo,angSep,projSep,massSep,MSgroup,SimGalID,HaloID,SimHaloMass,color\n");
      for(i = 1; i <= nsample; ++i){
        
        k = indx[i];
        igrp = group_member[k];
        j = group_index[igrp];

        if(k == j)
          theta = 0;
        if(k != j)
            theta = angular_separation(ra[k],dec[k],ra[j],dec[j]);

        fprintf(fp,"%d,%f,%f,%f,%f,%d,%f,%d,%f,%f,%f,%f,%f,%f,%d,%d,%f,%d\n", k, ra[k], dec[k], redshift[k]/speedOfLight, m_stellar[k], igrp, prob_total[k], group_center[igrp], mass[group_center[igrp]], rad[group_center[igrp]],theta, theta/angRad[group_center[igrp]], distToBig[igrp],group_mass[igrp],GalID[k],HaloID[k],m_halo[k],color_flag[k]);

      }
      fclose(fp);
      endout = get_msec() - startout;
      zoneout = (float)endout / 1000.;
    }
  }

  endall = get_msec() - startall;
  printf("%.3f sec\n", (float)endall/1000.); 

  for(i = 1; i <= niter_max; ++i){
    zone1t += zone1[i];
    zone2t += zone2[i];
    zone3t += zone3[i];
  }

  zone1t /= niter_max;
  zone2t /= niter_max;
  zone3t /= niter_max;

  printf("Initialization: %f seconds.\nIteration reset %f seconds (avg).\nIterative groupfinding: %f seconds (avg).\nSorting + sham: %f seconds (avg).\nOutput %f seconds.\n",zoneinit, zone1t, zone2t, zone3t, zoneout);

  exit(0);
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

float angular_separation(float ra1, float dec1, float ra2, float dec2)
{
    // return atan((sqrt(cos(dec2)*cos(dec2)*sin(ra2-ra1)*sin(ra2-ra1) + 
     //    pow(cos(dec1)*sin(dec2) - sin(dec1)*cos(dec2)*cos(ra2-ra1),2.0)))/
    //     (sin(dec1)*sin(dec2) + cos(dec1)*cos(dec2)*cos(ra2-ra1)));
  return(acos(sin(dec1) * sin(dec2) + cos(dec1) * cos(dec2) * cos(ra1 - ra2)));
}

float radial_probability(float mass, float dr, float rad, float ang_rad){
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

// This function uses the kdtree library written by John Tsiombikas <nuclear@member.fsf.org>.

float find_satellites(int i, float *ra, float *dec, float *redshift, int *color_flag, float theta_max, float sigma, int *group_member, int *indx, int ngal, float radius, float mass, int igrp, float *m_stellar, float *nsat_cur, float *prob_total, void *kd, float cenDist, float range, float sw1, float sw2) {
  int j, k, setSize;
  float dx, dy, dz, theta, prob_ang, prob_rad, grpMass, p0, sat_red_weight;
  void *set;
  int *pch;
  double cen[3];
  double sat[3];
  
  // Set up.
  
  dy = 1;
  dx = 1;
  *nsat_cur = 0;
  grpMass = 0;

  // Use the k-d tree kd to identify the nearest galaxies to the central.
  
  
  //cenDist = distance_redshift(redshift[i]/speedOfLight);
  cen[0] = cenDist * cos(ra[i]) * cos(dec[i]);
  cen[1] = cenDist * sin(ra[i]) * cos(dec[i]); 
  cen[2] = cenDist * sin(dec[i]);

  // Nearest neighbour search should go out to about 4*sigma, the velocity dispersion of the SHAMed halo.
  
  //range = distance_redshift(((4*sigma))/speedOfLight);
  set = kd_nearest_range(kd, cen, range);
  setSize = kd_res_size(set);
  
//  chunk = CHUNKSIZE;

  // Set now contains the nearest neighbours within a distance range. Grab their info. 
  // Note that set will ALWAYS contain a node that is the same as the central galaxy (this is a quirk of the code--or me not using it properly--when computing nearest neighbour distances to a point that is already in the k-d tree). Make sure to reject this galaxy.
  
 //  #pragma omp parallel num_threads (24) shared(group_member, redshift,prob_total,m_stellar,grpMass, set, sat,i) private(j, k)
//    {
//    #pragma omp for schedule (dynamic, chunk)
   for(k = 1; k <= setSize; ++k){
    //while( !kd_res_end(set)) {

    // Get the data and position of the current result item. Data contains index of galaxy.

    pch = (int*)kd_res_item(set, sat);
    j = *pch;
    
    // Move to next item in set. (If this isn't done before the next check and the first dist value is 0, everything gets wonky.)

    kd_res_next(set);

    // Are we at the end of the kd-tree? If so, flip the flag.
    
     
    //  flag = 0;
    //  continue;
    //}

    // Skip if target galaxy is the same as the central (obviously).
    if(i == j){
      continue;
    }

    // Skip if already assigned to a central.
    if(group_member[j]){
      continue;
    }

    // dx = fabs(ra[i]-ra[j]);
  //      if(dx>4*theta_max)
  //        continue;
  //      dy = fabs(dec[i]-dec[j]);
  //      if(dy>4*theta_max)
  //        continue;
        dz = fabs(redshift[i] - redshift[j]);
        // if(dz>6*sigma)
        //  continue;
        //printf("%d %f %f %d %f %f\n",i,ra[i],dec[i],j,ra[j],dec[j]);
    theta = angular_separation(ra[i],dec[i],ra[j],dec[j]);
    if(theta > theta_max){
      continue;
    }

    // Now determine the probability of being a satellite (both projected onto the sky, and along the line of sight).
    
    prob_ang = radial_probability(mass, theta, radius, theta_max);
    prob_rad = exp(-dz*dz/(2*sigma*sigma))*speedOfLight/(rt2Pi*sigma);
    
    p0 = (1 - 1/(1 + prob_ang * prob_rad / 10));
    
    if(p0 < 0){
      printf("ZERO %e\n",p0);
      p0 = 0;
    }

    if(p0 > prob_total[j])
      prob_total[j] = p0;
    if(prob_total[j] > 1)
      prob_total[j] = 1;
    if(prob_ang*prob_rad < 10){
      continue;
    }
    
    // At this point the galaxy is a satellite. Assign the group ID number to it, and add its mass to the total group mass. Increase satellite counter by 1.

    group_member[j] = igrp;
    ;
    // Check the satellite galaxy's color_flag variable. If it is set to 1, then apply the red weight for satellites.
    if(color_flag[j] == 1)
      sat_red_weight = (log10(m_stellar[j]) * sw1) + sw2;
    else
      sat_red_weight = 1;
    
    grpMass += (m_stellar[j] * sat_red_weight);
    (*nsat_cur) += 1;
  }
//  }
   
  // Vmax correction.

  // dz = speedOfLight* fabs(redshift[i] - MINREDSHIFT);
  // vol_corr = 1-(0.5*erfc(dz/(root2*sigma)));
  // *nsat_cur /= vol_corr;
  // grpMass /= vol_corr;
  
  // dz = speedOfLight* fabs(redshift[i] - MAXREDSHIFT);
  // vol_corr = 1-(0.5*erfc(dz/(root2*sigma)));
  // *nsat_cur /= vol_corr;
  // grpMass /= vol_corr;
 
  return(grpMass);

}

int central_galaxy(int i, float *ra, float *dec, int *group_member, int ngal, int igrp, float radius, float *m_stellar, float nsat){

  int j, ii[100000], icnt, k, imax, l, icen;
  float pot[100000], theta, maxpot, angDist;
    
  // This function identifies the central galaxy in a group by identifying the most massive galaxy. 

  icen = i;

  // Begin with target group, and find all galaxies within it by scanning the group_member array.

  icnt = 0;
  for(i = 1; i <= ngal; ++i){
    if(group_member[i] == igrp){
      icnt++;
      ii[icnt] = i;
    }
  }

  // Perform a double sum of mass to get the 'potential' of each galaxy.

    for(k = 1; k <= icnt; ++k)
    pot[k] = 0;

  for(k = 1; k <= icnt; ++k){
      for(l = k+1; l <= icnt; ++l){
        i = ii[k];
      j = ii[l];
      angDist = angular_separation(ra[i], dec[i], ra[j], dec[j]);
      theta = sqrt((angDist * angDist) + ((radius * radius) / 400.0));
      pot[k] += (m_stellar[i] * m_stellar[j]) / theta;
      pot[l] += (m_stellar[i] * m_stellar[j]) / theta;
      }
      }

  // Which galaxy member has the highest potential? This is considered to be the central galaxy.

  maxpot = 0;
  for(k = 1; k <= icnt; ++k){
      i = ii[k];
        //printf("%d %d %d %f %f %f %f %d\n",k,i,icen,ra[i],dec[i],pot[k],maxpot,imax);
        if(pot[k] > maxpot){ 
          maxpot = pot[k]; 
          imax = i;
      }
    }

  return imax;

}

int iter_central_galaxy(int galID, float *ra, float *dec, float *redshift ,int *group_member, int ngal, int igrp, float *m_stellar, float *massDist){
  float xcen = 0., ycen = 0., zcen = 0., Msum = 0., *distToCoM, maxDist, cenDist, satDist;
  float maxMass;

  int goodGals[2], maxMassID;
  int i, k, l, *groupGals, icnt, IterCenID, maxDistID, *rejectGals;
  int maxDistLoc, remainingGals;
  
  float cen[3], sat[3];
  
    
  // The purpose of this function is to compute the iterative central galaxy for a group, in the style of Robotham et al. 2011.

  // This is done by iteratively computing the centre of mass for the group, and rejecting the galaxy furthest from it. Once there are only 2 remaining galaxies, the more massive galaxy is chosen to be the group centre.

  // Begin with target group, and find all galaxies within it by scanning the group_member array.

  icnt = 0;
  for(i = 1; i <= ngal; ++i){
    if(group_member[i] == igrp){
      icnt++;
      Msum += m_stellar[i];
    }
  }
  
  groupGals = ivector(1, icnt);

  icnt = 0;
  for(i = 1; i <= ngal; ++i){
    if(group_member[i] == igrp){
      icnt++;
      groupGals[icnt] = i;
    }
  }
  
  distToCoM = vector(1,icnt);
  rejectGals = ivector(1, icnt);
  for(i = 1; i <= icnt; ++i){
    rejectGals[i] = 0;
  }
  maxDistLoc = 1;
  remainingGals = icnt;

  // Start by computing the centre of mass for the full group.

  for(k = 1; k <= icnt; ++k){

    i = groupGals[k];
    
    cenDist = distance_redshift(redshift[i]/speedOfLight);
    cen[0] = cenDist * cos(ra[i]) * cos(dec[i]); // x
    cen[1] = cenDist * sin(ra[i]) * cos(dec[i]); // y
    cen[2] = cenDist * sin(dec[i]); // z
    
    xcen += m_stellar[i] * cen[0];
    ycen += m_stellar[i] * cen[1];
    zcen += m_stellar[i] * cen[2];

  }

  xcen /= Msum;
  ycen /= Msum;
  zcen /= Msum;

  // These are the xyz coords of the first centre of mass. Loop through group and throw away furthest galaxy. Do this until icnt drops to 2.

    while(remainingGals > 2){

      for(k = 1; k <= icnt; ++k){

        if(rejectGals[k] != 0)
        continue;

        i = groupGals[k];

        satDist = distance_redshift(redshift[i]/speedOfLight);
      sat[0] = satDist * cos(ra[i]) * cos(dec[i]); // x
      sat[1] = satDist * sin(ra[i]) * cos(dec[i]); // y
      sat[2] = satDist * sin(dec[i]); // z

      distToCoM[k] = sqrt((xcen - sat[0])*(xcen-sat[0]) + (ycen - sat[1])*(ycen-sat[1]) + (zcen - sat[2])*(zcen-sat[2]));

      }

      // Find maximum distance. Grab the ID and the location in the array where it is.

      maxDist = 0;
      for(i = 1; i <= icnt; ++i){
        k = groupGals[i];
        if(distToCoM[i] > maxDist){
          maxDist = distToCoM[i];
          maxDistID = k;
          maxDistLoc = i;
        }
      }

      // Discard furthest galaxy and decrease remainingGals by 1. Remove this galaxy's stellar mass from total stellar mass.

      --remainingGals;
      rejectGals[maxDistLoc] = maxDistID;
      Msum -= m_stellar[maxDistID];

      // Clean various bits and pieces for next step in loop.

      for(i = 1; i <= icnt; ++i){
        distToCoM[i] = 0.;
      }

      // Compute new centre of mass position.
    
      for(k = 1; k <= icnt; ++k){
        if(rejectGals[k] != 0)
          continue;

      i = groupGals[k];
      
      cenDist = distance_redshift(redshift[i]/speedOfLight);
      cen[0] = cenDist * cos(ra[i]) * cos(dec[i]); // x
      cen[1] = cenDist * sin(ra[i]) * cos(dec[i]); // y
      cen[2] = cenDist * sin(dec[i]); // z
      
      xcen += m_stellar[i] * cen[0];
      ycen += m_stellar[i] * cen[1];
      zcen += m_stellar[i] * cen[2];

    }

    xcen /= Msum;
    ycen /= Msum;
    zcen /= Msum;

      // printf("%d %d %d %d\n",icnt, remainingGals, counter-1, maxDistID);
  }

  // There should only be two galaxies remaining now. Pick the one with the greater stellar mass as the group center.

  l = 0;
  for(k = 1; k <= icnt; ++k){

    if(rejectGals[k] != 0)
        continue;

    i = groupGals[k];

    goodGals[l] = i;
    ++l;
  }

  if(m_stellar[goodGals[0]] > m_stellar[goodGals[1]]){
    IterCenID = goodGals[0];
  }
  if(m_stellar[goodGals[0]] < m_stellar[goodGals[1]]){
    IterCenID = goodGals[1];
  }
  
  // Finally, what is the distance between the most massive galaxy in the group and the new group center?

  maxMass = 0.;
  for(k = 1; k <= icnt; ++k){
    i = groupGals[k];
    if(m_stellar[i] > maxMass){
      maxMass = m_stellar[i];
      maxMassID = i;
    }
  }

  cenDist = distance_redshift(redshift[maxMassID]/speedOfLight);
  cen[0] = cenDist * cos(ra[maxMassID]) * cos(dec[maxMassID]); // x
  cen[1] = cenDist * sin(ra[maxMassID]) * cos(dec[maxMassID]); // y
  cen[2] = cenDist * sin(dec[maxMassID]); // z

  satDist = distance_redshift(redshift[IterCenID]/speedOfLight);
  sat[0] = satDist * cos(ra[IterCenID]) * cos(dec[IterCenID]); // x
  sat[1] = satDist * sin(ra[IterCenID]) * cos(dec[IterCenID]); // y
  sat[2] = satDist * sin(dec[IterCenID]); // z

  *massDist = sqrt((cen[0] - sat[0])*(cen[0]-sat[0]) + (cen[1] - sat[1])*(cen[1]-sat[1]) + (cen[2] - sat[2]) * (cen[2]-sat[2]));

  free_vector(distToCoM,1,icnt);
  free_ivector(groupGals,1,icnt);
  free_ivector(rejectGals,1,icnt-2);
  return IterCenID;

}


// static double dist_sq( double *a1, double *a2, int dims ) {
//   double dist_sq = 0, diff;
//   while( --dims >= 0 ) {
//     diff = (a1[dims] - a2[dims]);
//     dist_sq += diff*diff;
//   }
//   return dist_sq;
// }

float segvol(float lorad, float hirad, float lodec, float hidec, float lora, float hira){
  float h, area, vol;
  lodec *= pi/180;  
  hidec *= pi/180;
  h = sin(hidec) - sin(lodec);
  area = 2*pi*((hira-lora)/360)*h;
  vol = (((4./3.) * pi * pow(hirad, 3.)) * area/(4. * pi)) - (((4./3.) * pi * pow(lorad, 3.)) * area/(4. * pi));
  return(vol);
}

float segarea(float lorad, float hirad, float lodec, float hidec, float lora, float hira){
  float h, area;
  lodec *= pi/180;
  hidec *= pi/180;
  h = sin(hidec) - sin(lodec);
  area = 2*pi*((hira-lora)/360)*h;
  return(area);
}
