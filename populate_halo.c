#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "nrutil.h"

// This code populates halos from an external file (output from a simulation) with mock galaxies, assuming some HOD. Normally one would parametrize this, but this version will read in a HOD file.

int main(int argc, char **argv){

	FILE *fp;
	int haloID, i, j, k, l, n, count, nhalo;
	double *mass, posgal[3], velgal[3], *poshalo[3], *velhalo[3], *ncen, *nsat;
	float mMin;
	char *ff;
	char string[1000];

	// Read in the halo file.

	ff = "/home/users/ma5046/work/bolshoi_redux.dat";

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
	printf("%d\n",nhalo);

	// Load halo information into arrays / matrices.

	mass = vector(1, nhalo);
	ncen = vector(1, nhalo);
	nsat = vector(1, nhalo);
	poshalo = matrix(1,nhalo,1,3);
	velhalo = matrix(1,nhalo,1,3);

	i = 0;
	for(j = 2; j <= nhalo+1; ++j){
		++i;
		fscanf(fp, "%f %f %f %*d %f %f %f", &poshalo[1,i], &poshalo[2,i], &poshalo[3,i],&velhalo[1,i], &velhalo[2,i], &velhalo[3,i],&mass[i]);
	}	

	poshalo[1,i], poshalo[2,i], poshalo[3,i],velhalo[1,i], velhalo[2,i], velhalo[3,i],mass[i]

}