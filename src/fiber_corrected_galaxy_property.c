#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* find a galaxy that it NOT fiber collided that is within \pm 0.25 magnitudes of the galaxy in question,
 * and within \pm 0.05 in g-r color, and submit that galaxy's property.
 */
float fiber_corrected_galaxy_property(int indx, float *mag_r, float *mag_g, float *prop, int ngal, int *flag)
{
  int i,iter=0;
  float c1, c2;

  c2 = mag_g[indx] - mag_r[indx];
  if(c2>1.2)c2 = drand48()*0.2+0.45;
  //printf("HERE %d %f %f\n",indx,mag_r[indx],c2);
  //fflush(stdout);

  while(1)
    {
      iter++;
      if(iter>1E6){ printf("BREAKFLAG %f %f %f\n",c2,mag_g[indx],mag_r[indx]); break;}
      i = drand48()*ngal+1;
      if(flag[i])continue;
      if(fabs(mag_r[i]-mag_r[indx])>0.25)continue;
      c1 = mag_g[i] - mag_r[i];
      if(fabs(c1-c2)<0.05)break;
    }
  
  //printf("DONE %d\n",indx);
  //fflush(stdout);


  return prop[i];

}
