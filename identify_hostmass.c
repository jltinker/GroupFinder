#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

#define OMEGA_M 0.25
#define PI 3.141592741
#define RHO_CRIT 2.775E+11
#define DELTA_HALO 200
#define THIRD (1.0/3.0)
#define BOXSIZE 200.0

int main(int argc, char **argv)
{
  FILE *fp;
  char a,a1[10],a2[10],a3[10],aa[1000];
  int i, j, k,i1, j1, nmatched, idum, n1, n2, ngal, ngal_sample;
  int **galid1, **galid2, xdum;
  float **galpos1, **galpos2, *redshift2, *mass2, *emass2; 
  float *magr, *magg, *ba90, *sersic_n;
  float *d4k, *d4k_err;
  int **ipos1, **ipos2, flag, nrepeat=0;
  float r, rmin, volume, vol_corr;
  float *mag_r, *mag_g, dx, dy, dz, x1, theta, theta_max, rmax, *ra, *dec, 
    *redshift, *mass, *rad, *sigma, *ang_rad, *vmax, *nsat_indi, *ptemp, 
    *luminosity, *group_luminosity;
  int *icollided, *indx, *group_member,*group_index;

  float prob_rad, prob_ang, *prob_total, p, x2, x3, err, mass_tolerance = 0.001, prev_mass, xngrp;
  int imass, niter_max=10, igrp, ngrp, *group_center, *temp_group, ngrp_tmp, niter, nsat_tot; 
  double ndens_gal = 0;
  double nsat[160], nhalo[160], nsat_cur, nsub[160]; 
  double nsat2[160], nhalo2[160], mbar2[160], mbars[160]; 
  double mbar[160];

  float *H_delta, *Dn4k, *mhost, *x, *y, *z;
  float vz,sig,rproj;
  int *icen, jmin;
  long IDUM=555;
  float hostmass;


  //read in all the galaxies from the VAGC
  fp = openfile("/Users/tinker/cosmo/CENSAT_DECOMP/NEW_NBODY_MOCKS/environ_quench/subhalo_neig_200_0.10_19.5.txt");
  ngal = filesize(fp)-7;

  x = vector(1,ngal);
  y = vector(1,ngal);
  z = vector(1,ngal);
  mag_r = vector(1,ngal);
  mhost = vector(1,ngal);
  rad = vector(1,ngal);
  mass = vector(1,ngal);
  icen = ivector(1,ngal);

  for(i=1;i<=7;++i)
    fgets(aa,1000,fp);

  for(i=1;i<=ngal;++i) 
    {
      fscanf(fp,"%f %f %f %f %f %f %f %f %f %d %f %d %f %d",
	     &x[i],&y[i],&z[i],&x1,&x1,&vz,&mass[i],&mag_r[i],&icen[i],&x1,&x1,&idum,&mhost[i],&idum);
      rad[i] = pow(3*mass[i]/(4*PI*DELTA_HALO*OMEGA_M*RHO_CRIT),THIRD);
      fgets(aa,1000,fp);      
    }
  fclose(fp);
  fprintf(stderr,"Done reading in lss.dat %d\n",ngal);

  for(i=1;i<=ngal;++i)
    {
      if(icen[i])
	{
	  printf("0.0 0.0 0\n");
	  continue;
	}
      rmin = 100;
      for(j=0;j<=ngal;++j)
	{
	  if(icen[j]==0)continue;

	  dx = fabs(x[i]-x[j]);
	  if(dx>BOXSIZE) dx = BOXSIZE-dx;
	  if(dx<0) dx = dx+BOXSIZE;
	  if(dx>rad[j]*2.5)continue;

	  dy = fabs(y[i]-y[j]);
	  if(dy>BOXSIZE) dy = BOXSIZE-dy;
	  if(dy<0) dy = dy+BOXSIZE;
	  if(dy>rad[j]*2.5)continue;

	  dz = fabs(z[i]-z[j]);
	  if(dz>BOXSIZE) dz = BOXSIZE-dz;
	  if(dz<0) dz = dz+BOXSIZE;
	  if(dz>rad[j]*2.5)continue;

	  r = sqrt(dx*dx+dy*dy+dz*dz);
	  if(r<rmin){ rmin = r; jmin = j; rproj = sqrt(dx*dx + dy*dy);}
	}
      if(rmin>99)
	{
	  printf("0.0 0.0 0\n");
	  continue;
	}
      printf("%f %f %d %f\n",rproj/rad[jmin],rmin,jmin,rad[jmin]);
    } 

}
